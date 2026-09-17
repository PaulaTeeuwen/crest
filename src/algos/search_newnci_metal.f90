!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2023 Philipp Pracht
!
! crest is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! crest is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with crest.  If not, see <https://www.gnu.org/licenses/>.
!================================================================================!

subroutine crest_search_newnci_metal(env,tim)
!*******************************************************************
!* New runtype "nci_metal": WIP, based on search_newnci.f90's
!* iMTD-style NCI workflow, used here as a template.
!* The mainloop/metadynamics part now runs as-is (plain GFN-FF,
!* no UFF4MOF-specific bias yet) once the metal centers have been
!* identified, typed, and printed.
!*******************************************************************
  use crest_parameters,only:wp,stdout,aatoau
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use shake_module
  use iomod
  use utilities
  use cregen_interface
  use uff4mof_data
  use crest_cn_module,only:calculate_CN
  use miscdata,only:RCOV
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich,m
  logical :: pr,wr,doreturn
!===========================================================!
  type(calcdata) :: calc
  type(mddata) :: mddat
  type(shakedata) :: shk

  type(mddata),allocatable :: mddats(:)
  integer :: nsim,nallout

  real(wp) :: energy,gnorm
  real(wp),allocatable :: grad(:,:)
  character(len=:),allocatable :: ensnam
  integer :: nat,nall
  real(wp),allocatable :: eread(:)
  real(wp),allocatable :: xyz(:,:,:)
  integer,allocatable  :: at(:)
  logical :: dump,ex
  character(len=80) :: atmp,btmp,str
  logical :: multilevel(6)
  logical :: start,lower
!===========================================================!
  logical :: isTMetal    !> external function, see bondconstraint.f90
  logical :: isSBlockMetal !> external function, see bondconstraint.f90 (Li,Na,K,...,Be,Mg,Ca,...)
  integer :: nmetal
  character(len=2) :: elem
  character(len=8),allocatable :: uffopts(:)
  integer :: nuffopts
  logical :: alreadylisted
  character(len=8) :: uffchoice
  logical :: uffchosen
  real(wp),allocatable :: cnvals(:)
  real(wp),allocatable :: bondmat(:,:)
  integer :: icn
  character(len=8) :: autolabel,effective
  logical :: autofound
  logical,allocatable :: specmatch(:,:),atlist(:)
!>--- UFF4MOF-informed coordination-sphere restraints
  character(len=8),allocatable :: efflabel(:)
  logical,allocatable :: effrestrain(:)
  real(wp),allocatable :: effanglefc(:) !> per-atom angle-restraint force constant (<=0 = use default)
!>--- r1=/theta0=: fully user-supplied override, bypassing the UFF4MOF
!>    table lookup for that specific quantity (see parse_calcdata.f90)
  real(wp),allocatable :: effcustomr1(:),effcustomtheta0(:)
  real(wp) :: labelr1,labeltheta0
  logical :: labelfound,have_label,have_customr1,have_customtheta0
  real(wp) :: uffr1,ufftheta0,rcovj_raw,rij_uff,dnow
  logical :: uffr1found
  type(constraint) :: newcons
  integer :: nrestrained
!>--- theta0 (target coordination geometry) angle restraints
  integer,allocatable :: donorlist(:)
  integer :: ndonors,a,b
  real(wp) :: v1(3),v2(3),ang_now,ang_target
  integer :: nangrestrained
!>--- auto=true geometric tiebreak among CN-tied UFF4MOF candidates
  logical :: uffauto
!>--- force=true: allow an explicit uff4mof= label from a DIFFERENT
!>    element's UFF4MOF namespace (deliberate cross-element borrowing)
  logical :: uffforce
  real(wp) :: forcer1,forcetheta0
  logical :: forcefound
  character(len=8),allocatable :: tiedlabels(:)
  integer,allocatable :: gdonorlist(:)
  integer :: ntied,ngdonors,itied,bestidx
  integer :: nundecided
  real(wp) :: score,bestscore
  logical :: suggestion_ok
  character(len=8) :: suggestlabel
!===========================================================!
!>--- printout header
  write (stdout,*)
  write (stdout,'(10x,"┍",49("━"),"┑")')
  write (stdout,'(10x,"│",13x,a,13x,"│")') " CREST NCI-METAL SAMPLING  "
  write (stdout,'(10x,"┕",49("━"),"┙")')
  write (stdout,*)

  write (stdout,*) 'Runtype nci_metal selected.'
  write (stdout,*) '(WIP: workflow is being filled in incrementally)'

!===========================================================!
!>--- setup
  call env%ref%to(mol)
  write (stdout,*) 'Input structure:'
  call mol%append(stdout)
  write (stdout,*)

!>--- identify metal atoms in the structure
  nmetal = 0
  do i = 1,mol%nat
    if (isanymetal(mol%at(i))) nmetal = nmetal+1
  end do

!>--- geometric coordination numbers, used below to auto-detect a
!>    matching UFF4MOF type for each metal atom (see uff4mof.f90).
!>    The raw exponential CN (cnvals) is fractional and tends to pick
!>    up small contributions from second-shell atoms, so the discrete
!>    coordination number used for UFF4MOF matching is instead a
!>    thresholded bond count from the same routine's "bond" matrix.
  allocate (cnvals(mol%nat))
  call calculate_CN(mol%nat,mol%at,mol%xyz,cnvals,cntype='exp',bond=bondmat)

  write (stdout,'(1x,a,1x,i0)') 'Number of metal atoms found:',nmetal
  if (nmetal > 0) then
    write (stdout,'(1x,a)') 'Metal atom identity:'
    do i = 1,mol%nat
      if (isanymetal(mol%at(i))) then
        write (stdout,'(3x,a,i0,a,a2,a,i0,a)') &
        & 'atom ',i,': ',i2e(mol%at(i),'nc'),'  (CN = ',count(bondmat(1:mol%nat,i) > 0.5_wp),', from geometry)'
      end if
    end do
  end if
  write (stdout,*)

!>--- for each distinct metal element present, list the UFF4MOF atom
!>    types tabulated for it (see uff4mof.f90, from lammps_interface)
  if (nmetal > 0) then
    write (stdout,'(1x,a)') 'Available UFF4MOF atom types for the metals present:'
    do i = 1,mol%nat
      if (.not.isanymetal(mol%at(i))) cycle
      elem = i2e(mol%at(i),'nc')
!>--- skip elements already listed
      alreadylisted = .false.
      do j = 1,i-1
        if (isanymetal(mol%at(j))) then
          if (trim(i2e(mol%at(j),'nc')) == trim(elem)) alreadylisted = .true.
        end if
      end do
      if (alreadylisted) cycle

      call uff4mof_options(elem,uffopts,nuffopts)
      if (nuffopts > 0) then
        write (stdout,'(3x,a,a)',advance='no') trim(elem),' : '
        do k = 1,nuffopts
          write (stdout,'(a,1x)',advance='no') trim(uffopts(k))
        end do
        write (stdout,*)
        deallocate (uffopts)
      else
        write (stdout,'(3x,a,a)') trim(elem),' : (no UFF4MOF parameters tabulated)'
      end if
    end do
  end if

!>--- resolve every [[calculation.nci_metal]] atom spec (a single index,
!>    a list/range, or an element symbol) against the actual structure.
!>    Where specs overlap for the same atom, the LAST one added wins.
  allocate (specmatch(mol%nat,max(env%nci_metal%nassign,1)))
  specmatch = .false.
  do k = 1,env%nci_metal%nassign
    call get_atlist(mol%nat,atlist,trim(env%nci_metal%spec(k)),mol%at)
    specmatch(1:mol%nat,k) = atlist(1:mol%nat)
    if (.not.any(atlist)) then
      write (stdout,'(1x,a,a,a)') '**WARNING** [[calculation.nci_metal]] atom spec "', &
      & trim(env%nci_metal%spec(k)),'" matched no atoms in this structure; ignored'
    else
      do i = 1,mol%nat
        if (atlist(i).and..not.isanymetal(mol%at(i))) then
          write (stdout,'(1x,a,a,a,i0,a,a2,a)') &
          & '**WARNING** [[calculation.nci_metal]] atom spec "',trim(env%nci_metal%spec(k)), &
          & '" matched atom ',i,' (',i2e(mol%at(i),'nc'),'), which is not a metal atom; ignored'
        end if
      end do
    end if
    if (allocated(atlist)) deallocate (atlist)
  end do
  write (stdout,*)

!>--- per-atom UFF4MOF assignment: auto-detect from the geometric CN,
!>    then let an explicit [[calculation.nci_metal]] .toml entry (if
!>    any, and if valid) override that detection.
!>    Default is UNRESTRAINED (source=.false.), regardless of whether
!>    the type is ambiguous or resolves unambiguously on its own --
!>    restraint is always an explicit ask (via uff4mof=, auto=true, or
!>    restrain=true in a matching [[calculation.nci_metal]] entry),
!>    never a silent side effect of typing being possible. A matching
!>    spec still overrides this per-atom (see the assignment loop
!>    below), so uff4mof=/auto=true continue to imply restrain=true
!>    unless that same entry also explicitly sets restrain=false.
  allocate (efflabel(mol%nat))
  allocate (effrestrain(mol%nat),source=.false.)
  allocate (effanglefc(mol%nat),source=-1.0_wp)
  allocate (effcustomr1(mol%nat),source=-1.0_wp)
  allocate (effcustomtheta0(mol%nat),source=-1.0_wp)
  allocate (gdonorlist(mol%nat))
  nundecided = 0
  if (nmetal > 0) then
    write (stdout,'(1x,a)') 'UFF4MOF assignment (per metal atom):'
    do i = 1,mol%nat
      if (.not.isanymetal(mol%at(i))) cycle
      elem = i2e(mol%at(i),'nc')
      icn = count(bondmat(1:mol%nat,i) > 0.5_wp)
      call uff4mof_autodetect(elem,icn,autolabel,autofound)
!>--- last matching [[calculation.nci_metal]] spec (if any) wins, both
!>    for the UFF4MOF label itself and for whether to restrain it
      uffchosen = .false.
      uffauto = .false.
      uffforce = .false.
      do k = env%nci_metal%nassign,1,-1
        if (specmatch(i,k)) then
          uffchoice = env%nci_metal%uff(k)
          uffchosen = len_trim(uffchoice) > 0
          uffauto = env%nci_metal%auto(k)
          uffforce = env%nci_metal%force(k)
          effrestrain(i) = env%nci_metal%restrain(k)
          effanglefc(i) = env%nci_metal%angle_fc(k)
          effcustomr1(i) = env%nci_metal%r1(k)
          effcustomtheta0(i) = env%nci_metal%theta0(k)
          exit
        end if
      end do
      call uff4mof_options(elem,uffopts,nuffopts)

!>--- determine, independent of auto=true, whether this is a GENUINE
!>    CN-tie (>=2 tabulated types share this CN, e.g. Cu4+2/Cu3f2 at
!>    CN4) as opposed to simply nothing tabulated at all -- needed
!>    both for the auto=true tiebreak below and for the halt-if-
!>    undecided check after this loop. When there IS a genuine tie,
!>    also compute (unconditionally, not just under auto=true) which
!>    candidate the CURRENT geometry fits best, so that even a plain
!>    "CN-ambiguous" report -- and, if this ends up unresolved, the
!>    halt message after this loop -- can show a concrete suggestion.
      ntied = 0
      ngdonors = 0
      suggestion_ok = .false.
      if (.not.autofound) then
        call uff4mof_candidates_at_cn(elem,icn,tiedlabels,ntied)
        if (ntied >= 2) then
          do j = 1,mol%nat
            if (j == i.or.bondmat(i,j) < 0.5_wp) cycle
            select case (i2e(mol%at(j),'nc'))
            case ('N','O','S','P','F','Cl','Br','I')
              ngdonors = ngdonors+1
              gdonorlist(ngdonors) = j
            end select
          end do
          if (ngdonors >= 2) then
            bestidx = 1
            call uff4mof_get_entry(tiedlabels(1),uffr1,ufftheta0,uffr1found)
            bestscore = angle_fit_score(ufftheta0,gdonorlist,ngdonors,i)
            do itied = 2,ntied
              call uff4mof_get_entry(tiedlabels(itied),uffr1,ufftheta0,uffr1found)
              score = angle_fit_score(ufftheta0,gdonorlist,ngdonors,i)
              if (score < bestscore) then
                bestscore = score
                bestidx = itied
              end if
            end do
            suggestlabel = tiedlabels(bestidx)
            suggestion_ok = .true.
          end if
        end if
      end if

      write (stdout,'(3x,a,i0,a,a2,a,i0,a)',advance='no') &
      & 'atom ',i,' (',elem,'), CN=',icn,': '
      if (autofound) then
        write (stdout,'(a,a)',advance='no') 'auto-detected ',trim(autolabel)
      else if (uffauto.and..not.uffchosen) then
!>--- auto=true opt-in: actually USE the geometric suggestion computed
!>    above (comparing actual donor angles against each tied
!>    candidate's theta0) to resolve the tie.
        if (suggestion_ok) then
          autolabel = suggestlabel
          autofound = .true.
          write (stdout,'(a)',advance='no') 'auto=true, tied by CN ('
          do itied = 1,ntied
            if (itied > 1) write (stdout,'(a)',advance='no') '/'
            write (stdout,'(a)',advance='no') trim(tiedlabels(itied))
          end do
          write (stdout,'(a,a)',advance='no') '), geometry best matches ',trim(autolabel)
        else
          write (stdout,'(a)',advance='no') 'auto=true requested but no geometric tiebreak possible '// &
          & '(not enough donor atoms found around this metal)'
        end if
      else if (ntied >= 2) then
        write (stdout,'(a)',advance='no') 'CN-ambiguous, tied between ('
        do itied = 1,ntied
          if (itied > 1) write (stdout,'(a)',advance='no') '/'
          write (stdout,'(a)',advance='no') trim(tiedlabels(itied))
        end do
        write (stdout,'(a)',advance='no') ')'
        if (suggestion_ok) then
          write (stdout,'(a,a,a)',advance='no') ' -- geometry suggests ',trim(suggestlabel), &
          & ' (add auto=true to use this, or uff4mof="..." to pick manually)'
        end if
      else
        write (stdout,'(a)',advance='no') 'no unambiguous auto-detection'
      end if

      effective = ''
      if (uffchosen) then
        if (nuffopts > 0.and.any(uffopts(1:nuffopts) == uffchoice)) then
          write (stdout,'(a,a)',advance='no') ', user override ',trim(uffchoice)
          effective = uffchoice
        else if (uffforce) then
!>--- force=true: accept a label from a DIFFERENT element's namespace,
!>    as long as it exists SOMEWHERE in the table at all (catches
!>    typos even under force=true) -- deliberately borrows that other
!>    element's tabulated r1/theta0 for THIS atom. Always loudly
!>    flagged: this is real, deliberate physics substitution, not a
!>    normal auto-detected or same-element manual choice.
          call uff4mof_get_entry(uffchoice,forcer1,forcetheta0,forcefound)
          if (forcefound) then
            write (stdout,'(a,a,a,a,a)',advance='no') &
            & ', **WARNING** user override "',trim(uffchoice),'" is for a DIFFERENT element than ', &
            & trim(elem),' -- forced anyway (force=true): borrowing its r1/theta0 deliberately'
            effective = uffchoice
          else
            write (stdout,'(a,a,a)',advance='no') &
            & ', **WARNING** user override "',trim(uffchoice),'" is not a tabulated UFF4MOF label at all; ignored'
            effective = autolabel
          end if
        else
          write (stdout,'(a,a,a)',advance='no') &
          & ', **WARNING** user override "',trim(uffchoice),'" is not tabulated for this element; ignored'
          effective = autolabel
        end if
      else
        effective = autolabel
      end if

      if (len_trim(effective) > 0) then
        write (stdout,'(a,a)') ' -> effective: ',trim(effective)
      else
        write (stdout,'(a)') ' -> effective: none (please specify manually)'
      end if
      efflabel(i) = effective

!>--- r1=/theta0=: report a custom override up front, alongside the
!>    UFF4MOF label resolution (which may itself be "none" -- a custom
!>    override does not require any resolved label at all).
      if (effcustomr1(i) > 0.0_wp.or.effcustomtheta0(i) > 0.0_wp) then
        write (stdout,'(6x,a)',advance='no') 'custom override:'
        if (effcustomr1(i) > 0.0_wp) write (stdout,'(a,f6.3,a)',advance='no') ' r1=',effcustomr1(i),' A'
        if (effcustomtheta0(i) > 0.0_wp) write (stdout,'(a,f6.2,a)',advance='no') ' theta0=',effcustomtheta0(i),' deg'
        write (stdout,'(a)') ' (bypasses UFF4MOF table for this quantity)'
      end if

!>--- Halt ONLY if the user showed explicit intent to restrain THIS
!>    atom and that intent could not be honored (genuine CN-tie, still
!>    unresolved). effrestrain(i) defaults to .false. (see allocate
!>    above) and is only ever set .true. by a matching
!>    [[calculation.nci_metal]] entry that itself carried an explicit
!>    uff4mof=, auto=true, or restrain=true -- so effrestrain(i) being
!>    .true. here already IS the "explicit intent" signal, with no
!>    separate specmatch check needed. A metal with NO matching entry
!>    at all never halts and never gets restrained -- "setting
!>    nothing" means an ordinary, unrestrained nci_metal run, exactly
!>    as before this feature existed, regardless of whether its type
!>    happens to be ambiguous OR resolves unambiguously on its own.
!>    Only an explicit ask that goes unanswered halts. A custom r1/
!>    theta0 override also counts as resolved here, even with no
!>    UFF4MOF label at all -- the ambiguity being flagged is about
!>    which TABLE entry to use, and a custom override bypasses the
!>    table entirely.
      if (len_trim(effective) == 0.and.ntied >= 2.and.effrestrain(i) &
      &   .and.effcustomr1(i) <= 0.0_wp.and.effcustomtheta0(i) <= 0.0_wp) then
        nundecided = nundecided+1
      end if

      if (allocated(tiedlabels)) deallocate (tiedlabels)
      if (allocated(uffopts)) deallocate (uffopts)
    end do
  end if
  deallocate (gdonorlist)
  write (stdout,*)

!>--- Stop here, before any MTD is run, if one or more metal atoms
!>    have a genuinely ambiguous UFF4MOF type (>=2 tabulated
!>    candidates share the observed CN) that was never resolved.
!>    Proceeding anyway would spend the (often hours-long) MTD budget
!>    exploring a coordination geometry that was never actually
!>    decided -- see the per-atom "CN-ambiguous" printout above for
!>    which candidates were tied. To proceed:
!>      - give an explicit  uff4mof = "..."  label, or
!>      - set  auto = true  (works if the metal has >=2 identifiable
!>        donor atoms already in range, so the geometric tiebreak has
!>        something to compare against), or
!>      - set  restrain = false  to explicitly acknowledge the
!>        ambiguity and proceed without deciding/restraining it.
  if (nundecided > 0) then
    write (stdout,'(1x,a,i0,a)') '**STOPPING**: ',nundecided, &
    & ' metal atom(s) above have a genuinely CN-ambiguous UFF4MOF type that was not resolved.'
    write (stdout,'(1x,a)') 'Add an explicit uff4mof="...", set auto=true, or set restrain=false '// &
    & 'in a [[calculation.nci_metal]] entry for each flagged atom, then rerun.'
    return
  end if

!>--- UFF4MOF-informed coordination-sphere detection and restraints.
!>    crest's own cn_to_bond/calculate_CN (bondmat above) uses one
!>    generic, inflated (x4/3) covalent radius per element -- fine for
!>    an ordinary bond, but its 50%-bonded cutoff sits well *inside* the
!>    range where a real, if weak or long, dative bond can still live,
!>    so a hemilabile or secondary donor can end up "not bonded" and
!>    then, since GFN-FF's topology is static (fixed once at setup, see
!>    gfnff_setup.f90), simply drift away over the course of an MTD.
!>    Here we instead use the *specific* UFF4MOF type resolved above
!>    (r1, per oxidation state/CN, from uff4mof.f90) to estimate a
!>    metal-ligand equilibrium length, search a generous window around
!>    it for candidate donor atoms, and -- unless the user disabled it
!>    via restrain=false -- add a soft flat-bottom distance restraint
!>    (bondrange: free to move near equilibrium, walled off before full
!>    dissociation) so the assigned coordination sphere survives the
!>    MTD bias potential regardless of what GFN-FF's own topology saw.
!>    This only ever *adds* restraints on top of GFN-FF; it does not
!>    touch GFN-FF's own topology or try to force an exact coordination
!>    number -- only atoms that are geometrically plausible donors (see
!>    donor-element list and the search window below) are restrained.
!>    Distance is handled by a soft bondrange restraint (free to move
!>    near equilibrium, walled off before full dissociation). Target
!>    *geometry* is handled separately below: UFF4MOF's theta0 is the
!>    idealized angle for the metal's assigned type (e.g. 90 deg for
!>    both square-planar Pd4+2 and octahedral Zn6+2), so every donor
!>    PAIR is restrained toward whichever of {theta0,180 deg} its
!>    CURRENT angle is closer to -- i.e. the starting geometry decides
!>    which pairs are "cis" (-> theta0) vs "trans" (-> 180 deg), and
!>    that assignment is then held for the rest of the run. This is
!>    what actually lets e.g. Pd4+2 (square-planar) be forced against
!>    GFN-FF's own electronic preference for tetrahedral, since GFN-FF
!>    (a classical force field) has no d-electron-count-driven geometry
!>    preference of its own to fall back on.
  if (nmetal > 0) then
    nrestrained = 0
    nangrestrained = 0
    allocate (donorlist(mol%nat))
    write (stdout,'(1x,a)') 'Coordination-sphere restraints (from UFF4MOF r1/theta0):'
    do i = 1,mol%nat
      if (.not.isanymetal(mol%at(i))) cycle
      have_label = len_trim(efflabel(i)) > 0
      have_customr1 = effcustomr1(i) > 0.0_wp
      have_customtheta0 = effcustomtheta0(i) > 0.0_wp
!>--- a metal-ligand bond-radius estimate (from the resolved UFF4MOF
!>    label's r1, or a user-supplied custom r1) is required just to
!>    search for donor atoms at all -- without one, there is nothing to
!>    restrain against, custom theta0 override or not (see r1/theta0
!>    docs in parse_calcdata.f90).
      if (.not.have_label.and..not.have_customr1) cycle
      labelfound = .false.
      if (have_label) call uff4mof_get_entry(efflabel(i),labelr1,labeltheta0,labelfound)
      if (have_customr1) then
        uffr1 = effcustomr1(i)
      else if (labelfound) then
        uffr1 = labelr1
      else
        cycle !> label resolved to an unusable entry and no custom r1 given
      end if
      if (have_customtheta0) then
        ufftheta0 = effcustomtheta0(i)
      else if (labelfound) then
        ufftheta0 = labeltheta0
      else
        ufftheta0 = 0.0_wp !> no angle target available -- distance restraint only (gated below)
        if (have_customr1) write (stdout,'(3x,a,i0,a)') &
        & '**NOTE** atom ',i,': custom r1 given but no theta0 (from override or UFF4MOF label) -- '// &
        & 'distance restraint only, no angle restraint'
      end if
      ndonors = 0
      do j = 1,mol%nat
        if (j == i) cycle
        if (isanymetal(mol%at(j))) cycle !> metal-metal contacts not handled here
!>--- restrict candidate donors to the usual Lewis-basic heteroatoms;
!>    a nearby C or H is essentially never itself the coordinating atom
        select case (i2e(mol%at(j),'nc')) !> i2e(...,'nc') is title-case, e.g. "Cl","Br"
        case ('N','O','S','P','F','Cl','Br','I')
        case default
          cycle
        end select
        rcovj_raw = RCOV(mol%at(j))*0.75_wp !> undo crest's own x4/3 CN-counting inflation
        rij_uff = uffr1+rcovj_raw/aatoau    !> both terms in Angstrom
        dnow = sqrt(sum((mol%xyz(:,i)-mol%xyz(:,j))**2))/aatoau !> Angstrom
!>--- generous detection window: catches weak/long dative contacts
!>    cn_to_bond's fixed cutoff would miss, without grabbing atoms that
!>    are simply nearby for unrelated geometric reasons
        if (dnow > 1.4_wp*rij_uff) cycle
        if (bondmat(i,j) < 0.5_wp) then
          write (stdout,'(3x,a,i0,a,i0,a,a2,a,f5.2,a,f5.2,a)') &
          & 'atom ',i,'-',j,' (',i2e(mol%at(j),'nc'),'): r=',dnow, &
          & ' A, not in cn_to_bond''s bondmat but within UFF4MOF window (target ',rij_uff,' A)'
        end if
        if (.not.effrestrain(i)) cycle
        call newcons%bondrangeconstraint(i,j,1.25_wp*rij_uff*aatoau,0.0_wp)
        call env%calc%add(newcons)
        nrestrained = nrestrained+1
        ndonors = ndonors+1
        donorlist(ndonors) = j
      end do
!>--- angle restraints among this metal's donor set. Two regimes:
!>    (1) theta0 ~ 90 deg (square-planar CN4 or octahedral CN6 -- the
!>        only UFF4MOF geometries with a genuine ANTIPODAL/trans
!>        relationship): each pair is assigned cis(theta0)/trans(180)
!>        role from the STARTING geometry, whichever it's closer to,
!>        then held there. Targeting 180 deg for SOME pairs here is not
!>        arbitrary -- it is a GEOMETRIC NECESSITY, not a choice: 4
!>        ligands genuinely arranged in a plane (or octahedrally)
!>        cannot have all 6 pairwise donor-metal-donor angles be 90 deg
!>        simultaneously -- the two ligands diagonally across from each
!>        other are unavoidably 180 deg apart. This is the mirror image
!>        of the tetrahedral case in regime (2) below, which has NO
!>        such geometric requirement at all.
!>        THIS IS ALSO WHY A RESTRAINED RUN CANNOT ISOMERIZE (cis <->
!>        trans): the role assignment is per SPECIFIC ATOM PAIR (e.g.
!>        "Cl2-Cu-N4 gets the 180 target"), decided once from the input
!>        and then held there by a continuous harmonic pull for the
!>        whole run. Both isomers have the same pattern of angles (4
!>        pairs at ~90, 2 at ~180) -- what differs between cis and
!>        trans is WHICH SPECIFIC PAIRS hold which role. Isomerizing
!>        would require reassigning that per-pair mapping mid-run,
!>        which never happens -- so cis started here can only ever stay
!>        cis, and vice versa (see cis-CuComplex/, trans-CuComplex/,
!>        cisplatin/, transplatin/ in final_tests/, and their README).
!>    (2) any other theta0 (e.g. 109.47 deg tetrahedral, 120 deg
!>        trigonal-planar): these geometries have NO 180-deg pairs at
!>        all -- ALL donor pairs are equivalent, so EVERY pair targets
!>        theta0 directly, no cis/trans split. Getting this wrong is
!>        not a minor issue: splitting by "closer to 180 or theta0" for
!>        a tetrahedral target permanently locks whichever pairs
!>        started near 180 (inherited from a square-planar input) to
!>        STAY at 180 forever -- a real tetrahedral center has no such
!>        pairs, so that lock is geometrically IMPOSSIBLE to satisfy
!>        together with the other pairs' theta0 target, regardless of
!>        angle_fc magnitude. Confirmed empirically: even a 200x
!>        stronger-than-default force constant (angle_fc=2.0, costing
!>        several hundred kcal/mol of unrelieved strain) could not
!>        relocate a square-planar-started Cu(II) structure to
!>        tetrahedral -- the final structure always had exactly 2 of
!>        6 pairs immovably stuck at 180.0 deg, a target that cannot
!>        coexist with 4 other pairs at 109.47 deg in real 3D space.
!>
!>    TODO(user-requested, not yet implemented): even within regime
!>    (1), the cis/trans ROLE of each pair is only ever inferred from
!>    the input geometry, so whatever isomer the input happens to be
!>    built as is what gets locked in and held for the rest of the
!>    run -- there is currently no way to ask for the opposite isomer
!>    (e.g. force trans when the input is cis). This matters
!>    chemically: trans is sometimes the genuinely lower-energy isomer
!>    (cf. cisplatin vs. transplatin), so geometry-inferred roles alone
!>    can't explore/confirm that. A later pass should let the user
!>    override which donor pairs are cis vs. trans (e.g. via
!>    [[calculation.nci_metal]]) instead of only ever trusting the
!>    starting structure.
!>
!>    CAVEAT for regime (2), discovered testing r1=/theta0= custom
!>    overrides (see parse_calcdata.f90): a UNIFORM target applied to
!>    ALL C(ndonors,2) pairs is only geometrically self-consistent for
!>    a handful of special values -- for exactly 4 donors, 109.47 deg
!>    (the regular tetrahedron) is the UNIQUE value where all 6
!>    pairwise angles can be simultaneously satisfied; any Td symmetry-
!>    breaking distortion that changes even one pair's angle away from
!>    109.47 necessarily changes at least one OTHER pair's angle too
!>    (this is exactly why regime (1) needs a genuine cis/trans SPLIT
!>    instead of one uniform target: square-planar/octahedral are the
!>    n=4/n=6 examples where the self-consistent solution is two
!>    distinct values, not one). Picking an "in-between" uniform theta0
!>    for 4 donors (e.g. 95 deg, intended as a distorted-tetrahedral
!>    target) asks for something no 4-donor geometry can actually
!>    satisfy across all 6 pairs at once -- confirmed empirically: at
!>    theta0=95 the restrained structure stayed within ~0.5 deg of
!>    109.47 regardless of angle_fc (tested 0.01 through 0.3; higher
!>    values caused CREGEN topology-check crashes on this system,
!>    unrelated to the restraint itself), i.e. the restraint pulled
!>    negligibly off GFN-FF's own self-consistent tetrahedral optimum.
!>    A genuinely distorted (non-uniform) 4-coordinate target needs a
!>    per-pair role split analogous to regime (1)'s cis/trans logic,
!>    which does not currently exist for non-90-deg centers -- until it
!>    does, theta0= for a 4-donor center should be set to a value that
!>    is ACTUALLY achievable uniformly (109.47 for tetrahedral being
!>    the only nontrivial one), not an arbitrary "somewhere between
!>    square-planar and tetrahedral" compromise value.
      if (ndonors >= 2.and.ufftheta0 > 0.0_wp) then
        do a = 1,ndonors-1
          do b = a+1,ndonors
            v1 = mol%xyz(:,donorlist(a))-mol%xyz(:,i)
            v2 = mol%xyz(:,donorlist(b))-mol%xyz(:,i)
            ang_now = acos(max(-1.0_wp,min(1.0_wp,dot_product(v1,v2)/(norm2(v1)*norm2(v2)))))*(180.0_wp/pi)
            if (abs(ufftheta0-90.0_wp) < 1.0_wp) then
!>--- regime (1): square-planar/octahedral -- genuine cis/trans roles exist
              if (abs(ang_now-180.0_wp) < abs(ang_now-ufftheta0)) then
                ang_target = 180.0_wp
              else
                ang_target = ufftheta0
              end if
            else
!>--- regime (2): tetrahedral/trigonal/etc -- no trans role exists,
!>    every pair targets theta0 directly
              ang_target = ufftheta0
            end if
            if (effanglefc(i) > 0.0_wp) then
              call newcons%angleconstraint(donorlist(a),i,donorlist(b),ang_target,k=effanglefc(i))
            else
              call newcons%angleconstraint(donorlist(a),i,donorlist(b),ang_target)
            end if
            call env%calc%add(newcons)
            nangrestrained = nangrestrained+1
          end do
        end do
      end if
    end do
    deallocate (donorlist)
    if (nrestrained > 0) then
      write (stdout,'(3x,i0,a)') nrestrained,' metal-ligand distance restraint(s) added'
      write (stdout,'(3x,i0,a)') nangrestrained,' donor-metal-donor angle restraint(s) added (theta0/trans)'
    else
      write (stdout,'(3x,a)') 'none added (either no resolved UFF4MOF type, or restrain=false)'
    end if
  end if
  write (stdout,*)

!>--- saftey terminations
  call crest_sampling_skip(env,doreturn)
  if (doreturn) return

!>--- sets the MD length according to a flexibility measure
  call md_length_setup(env)
!>--- create the MD calculator saved to env
  call env_to_mddat(env)

!===========================================================!
  if (env%performMTD) then
!>--- (optional) calculate a short 1ps test MTD to check settings
    call tim%start(1,'Trial metadynamics (MTD)')
    call trialmd(env)
    call tim%stop(1)
    if (env%iostatus_meta .ne. 0) return
  end if

!===========================================================!
!>--- Start mainloop
  env%nreset = 0
  start = .true.
  MAINLOOP: do
    call printiter
    if (.not.start) then
!>--- clean Dir for new iterations, but leave iteration backup files
      call clean_V2i
      env%nreset = env%nreset+1
    else
!>--- at the beginning, wipe directory clean
      call V2cleanup(.false.)
    end if
!===========================================================!
!>--- Meta-dynamics loop
    mtdloop: do i = 1,env%Maxrestart

      write (stdout,*)
      write (stdout,'(1x,a)') '------------------------------'
      write (stdout,'(1x,a,i0)') 'Meta-Dynamics Iteration ',i
      write (stdout,'(1x,a)') '------------------------------'

      nsim = -1 !>--- enambles automatic MTD setup in init routines
      call crest_search_multimd_init(env,mol,mddat,nsim)
      allocate (mddats(nsim),source=mddat)
      call crest_search_multimd_init2(env,mddats,nsim)

      call tim%start(2,'Metadynamics (MTD)')
      call crest_search_multimd(env,mol,mddats,nsim)
      call tim%stop(2)
!>--- a file called crest_dynamics.trj.xyz should have been written
      ensnam = 'crest_dynamics.trj.xyz'
!>--- deallocate for next iteration
      if (allocated(mddats)) deallocate (mddats)

!==========================================================!
!>--- Reoptimization of trajectories
      call tim%start(3,'Geometry optimization')
      call optlev_to_multilev(env%optlev,multilevel)
      call crest_multilevel_oloop(env,ensnam,multilevel,i)
      call tim%stop(3)
      if (env%iostatus_meta .ne. 0) return

!>--- save the CRE under a backup name
      call checkname_xyz(crefile,atmp,str)
      call checkname_xyz('.cre',str,btmp)
      call rename(atmp,btmp)
!>--- save cregen output
      call checkname_tmp('cregen',atmp,btmp)
      call rename('cregen.out.tmp',btmp)

!=========================================================!
!>--- cleanup after first iteration and prepare next
      if (i .eq. 1.and.start) then
        start = .false.
!>-- obtain a first lowest energy as reference
        env%eprivious = env%elowest
!>-- remove the two extreme-value MTDs
        if (.not.env%readbias.and.env%runver .ne. 33.and. &
        &   env%runver .ne. 787878) then
          env%nmetadyn = env%nmetadyn-2
        end if
!>-- the cleanup
        call clean_V2i
!>-- and always do two cycles of MTDs
        cycle mtdloop
      end if
!=========================================================!
!>--- Check for lowest energy
      call elowcheck(lower,env)
      if (.not.lower) then
        exit mtdloop
      end if
    end do mtdloop
!=========================================================!
!>--- collect all ensembles from mtdloop and merge
    write (stdout,*)
    write (stdout,'(''========================================'')')
    write (stdout,'(''           MTD Simulations done         '')')
    write (stdout,'(''========================================'')')
    write (stdout,'(1x,''Collecting ensmbles.'')')
!>-- collecting all ensembles saved as ".cre_*.xyz"
    call collectcre(env)
    call newcregen(env,0)
    call checkname_xyz(crefile,atmp,btmp)
!>--- remaining number of structures
    call remaining_in(atmp,env%ewin,nallout)

!==========================================================!
!>--- exit mainloop
    exit MAINLOOP
  end do MAINLOOP

!==========================================================!
!>--- final ensemble optimization
  write (stdout,'(/)')
  write (stdout,'(3x,''================================================'')')
  write (stdout,'(3x,''|           Final Geometry Optimization        |'')')
  write (stdout,'(3x,''================================================'')')
  call tim%start(3,'Geometry optimization')
  call checkname_xyz(crefile,atmp,str)
  call crest_multilevel_wrap(env,trim(atmp),0)
  call tim%stop(3)
  if (env%iostatus_meta .ne. 0) return

!==========================================================!
!>--- print CREGEN results and clean up Directory a bit
  write (stdout,'(/)')
  call smallhead('Final Ensemble Information')
  call V2terminating()

!========================================================================================!
contains
!========================================================================================!
!> nci_metal's own notion of "metal": transition metals (isTMetal, the
!> original d-block-only definition from bondconstraint.f90) plus the
!> alkali/alkaline-earth s-block metals (isSBlockMetal) -- both form
!> real, UFF4MOF-tabulated coordination environments. Kept local to this
!> routine rather than widening isTMetal itself, see isSBlockMetal's docs.
!========================================================================================!
  logical function isanymetal(z)
    implicit none
    integer,intent(in) :: z
    isanymetal = isTMetal(z).or.isSBlockMetal(z)
  end function isanymetal

!========================================================================================!
!> Sum of squared deviations between the ACTUAL donor-metal-donor angles
!> (donors(1:ndon), around atom `center`, from the current mol%xyz) and
!> whichever of {theta0test,180deg} each angle is closer to. Lower is a
!> better geometric fit to a candidate UFF4MOF type's theta0 -- used to
!> break a CN-tie between e.g. Cu4+2 (theta0=90) and Cu3f2 (theta0=109.47)
!> when the user opts in via auto=true. See the assignment loop above.
!========================================================================================!
  real(wp) function angle_fit_score(theta0test,donors,ndon,center)
    implicit none
    real(wp),intent(in) :: theta0test
    integer,intent(in) :: donors(:),ndon,center
    integer :: aa,bb
    real(wp) :: vv1(3),vv2(3),angnow,angtarget,dev
    angle_fit_score = 0.0_wp
    do aa = 1,ndon-1
      do bb = aa+1,ndon
        vv1 = mol%xyz(:,donors(aa))-mol%xyz(:,center)
        vv2 = mol%xyz(:,donors(bb))-mol%xyz(:,center)
        angnow = acos(max(-1.0_wp,min(1.0_wp,dot_product(vv1,vv2)/(norm2(vv1)*norm2(vv2)))))*(180.0_wp/pi)
        if (abs(angnow-180.0_wp) < abs(angnow-theta0test)) then
          angtarget = 180.0_wp
        else
          angtarget = theta0test
        end if
        dev = angnow-angtarget
        angle_fit_score = angle_fit_score+dev*dev
      end do
    end do
  end function angle_fit_score
!========================================================================================!
end subroutine crest_search_newnci_metal
!========================================================================================!
