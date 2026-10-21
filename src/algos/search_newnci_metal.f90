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
  use irmsd_module,only:irmsd,rmsd_cache
  use dynamics_module
  use shake_module
  use iomod
  use utilities
  use cregen_interface
  use uff4mof_data
  use crest_cn_module,only:calculate_CN
  use miscdata,only:RCOV
  use isomer_generate_mod,only:generate_isomer_structures
  use optimize_module,only:optimize_geometry
  use crest_restartlog
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
!>--- restart checkpoint (crest.restart) -- unseeded path only, see
!>    the detection block below for why the seeded path is excluded
  type(restart_data) :: rdat
  logical :: do_restart,skip_collect,fex
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
  integer :: total_cn    !> summed CN over all metal atoms, used to scale
                          !> MTD length by ligand/donor count (see below)
  real(wp) :: ligand_scale
  real(wp) :: ligflex
  integer,allocatable :: allinclude(:)
  character(len=8) :: autolabel,effective
  logical :: autofound
  logical,allocatable :: specmatch(:,:),atlist(:)
!>--- UFF4MOF-informed coordination-sphere restraints
  character(len=8),allocatable :: efflabel(:)
  logical,allocatable :: effrestrain(:)
  logical,allocatable :: effseed(:) !> per-atom resolved "seed" flag, see nci_metal_settings.f90
!>--- isomer seeding ("seed=true"): enumerate every symmetry-distinct
!>    arrangement of this metal's own ligands and run the full existing
!>    Kpush/alpha grid on EACH one (isomer_generate_mod), not just the
!>    input's own arrangement
  logical :: do_seed
  integer :: seed_metal_idx,n_seed_structs,n_total_structs,iblock
  real(wp),allocatable :: seed_structures(:,:,:)
  character(len=32) :: seed_geomname
  logical :: seed_ok
  type(calcdata),allocatable :: seedcalcs(:) !> one independently-initialized
    !> calculator per seed block (1=original input, 2..=each generated
    !> isomer), so GFN-FF topology is perceived separately per block
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
  logical :: is_carbon_donor,skip_distance,have_metal_r1
  integer :: ncarbondonors,nnodistdonors
!>--- theta0 (target coordination geometry) angle restraints
  integer,allocatable :: donorlist(:)
  integer :: ndonors,a,b
  real(wp) :: ang_target
  integer :: nangrestrained
!>--- flat-bottom half-width for the donor-metal-donor angle restraints
!>    below (theta0 and 180-deg trans targets alike): small enough to stay
!>    well inside any real gap between competing UFF4MOF coordination
!>    geometries (e.g. 90 vs. 109.47 deg), so it never leaves GFN-FF room
!>    to relax into a competing geometry unopposed, but wide enough to
!>    absorb ordinary bending-vibration noise instead of pulling on every
!>    single sampled/optimized structure regardless of how close it
!>    already is to the target.
  real(wp),parameter :: angletol = 4.0_wp !> degrees
!>--- regime (1) (theta0~90 deg) cis/trans role assignment: which
!>    SPECIFIC donor pairs get the antipodal 180-deg target, chosen as
!>    the best-fitting perfect matching over the CURRENT geometry (see
!>    best_trans_matching below) rather than a naive per-pair greedy
!>    threshold -- the greedy version degenerates for a starting
!>    geometry with no directional bias (e.g. exactly tetrahedral; see
!>    tetraplatin/ in final_tests/ and its README entry).
  integer,allocatable :: transpartner(:)
!>--- regime (3) (ndonors==5, theta0~120 deg) axial-pair role
!>    assignment: which SPECIFIC donor pair gets the 180-deg axial-axial
!>    target, chosen as the best-fitting of all C(5,2)=10 candidate
!>    pairs over the CURRENT geometry (see best_axial_pair_tbp below).
!>    Local indices into donorlist(:), same convention as transpartner.
  integer :: tbp_axial_a,tbp_axial_b
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

!>--- resolve distinguish_enantiomers HERE, not at TOML-parse time: it's
!>    a plain root-level toml key that may be parsed before OR after
!>    "runtype" depending on file order -- resolving it here, after all
!>    parsing has finished, is immune to that ordering. Default (false)
!>    leaves iinversion at CREST's own normal default (0, auto merge);
!>    opting in sets 2 (distinguish), so e.g. a Bailar-twist Delta/Lambda
!>    pair survives as two separate ensemble entries instead of merging.
  if (env%distinguish_enantiomers) then
    env%iinversion = 2
    write (stdout,'(1x,a)') 'distinguish_enantiomers = true: metal-centered '// &
    & 'enantiomers (e.g. Delta/Lambda) will be kept as distinct ensemble entries, '// &
    & 'not merged via iRMSD.'
  end if

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
!>--- seed defaults to false, same convention as restrain above: an
!>    opt-in feature should not silently run (and silently exercise
!>    isomer-generation + the enantiomer-duplicate filter) on every
!>    nci_metal call just because a metal atom was detected.
  allocate (effseed(mol%nat),source=.false.)
  allocate (effanglefc(mol%nat),source=-1.0_wp)
  allocate (effcustomr1(mol%nat),source=-1.0_wp)
  allocate (effcustomtheta0(mol%nat),source=-1.0_wp)
  allocate (gdonorlist(mol%nat))
  nundecided = 0
  total_cn = 0
  if (nmetal > 0) then
    write (stdout,'(1x,a)') 'UFF4MOF assignment (per metal atom):'
    do i = 1,mol%nat
      if (.not.isanymetal(mol%at(i))) cycle
      elem = i2e(mol%at(i),'nc')
      icn = count(bondmat(1:mol%nat,i) > 0.5_wp)
      total_cn = total_cn+icn
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
          effseed(i) = env%nci_metal%seed(k)
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
            call print_uff_typo_hint(elem,uffopts,nuffopts)
            effective = autolabel
          end if
        else
!>--- label belongs to a DIFFERENT element's UFF4MOF namespace and
!>    force=true was NOT given -- this is the accidental-typo/wrong-
!>    element-guess case the halt-before-MTD safety net (see below)
!>    does not otherwise catch, since a label WAS given (so this atom
!>    never counts as "undecided"). Reject it loudly rather than
!>    silently falling back to auto-detection, and point at both
!>    escape routes: force=true (if the cross-element borrow was
!>    actually intentional) and the list of labels genuinely tabulated
!>    for THIS element (if it was a typo/mistaken guess).
          write (stdout,'(a,a,a,a,a)',advance='no') &
          & ', **WARNING** user override "',trim(uffchoice),'" is not a tabulated UFF4MOF type for ', &
          & trim(elem),'; ignored (set force=true to deliberately borrow it from another element anyway)'
          call print_uff_typo_hint(elem,uffopts,nuffopts)
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
      labelfound = .false.
      if (have_label) call uff4mof_get_entry(efflabel(i),labelr1,labeltheta0,labelfound)
!>--- a metal-ligand bond-radius estimate (UFF4MOF label's r1, or a
!>    user-supplied custom r1) is needed for a DISTANCE restraint and
!>    for the donor-detection WINDOW search below -- but NOT for an
!>    angle restraint, which is pure direction, no bond length
!>    involved. have_metal_r1=false (no label, no custom r1 -- e.g. no
!>    tabulated UFF4MOF entry for this CN at all, same situation as the
!>    carbon-donor case below) no longer skips this metal outright as
!>    long as a theta0 (custom, or via the "geometry=" name shorthand,
!>    which resolves to the same effcustomtheta0) is available: donor
!>    detection falls back to bondmat (same mechanism as carbon
!>    donors), and every donor gets angle-restraint geometry awareness
!>    with no distance restraint, instead of getting nothing at all.
      have_metal_r1 = have_label.or.have_customr1
      if (.not.have_metal_r1.and..not.have_customtheta0) cycle
      if (have_metal_r1) then
        if (have_customr1) then
          uffr1 = effcustomr1(i)
        else if (labelfound) then
          uffr1 = labelr1
        else
          have_metal_r1 = .false. !> label resolved to an unusable entry
        end if
      end if
      if (.not.have_metal_r1.and..not.have_customtheta0) cycle !> re-check after the labelfound fallback above
      if (have_customtheta0) then
        ufftheta0 = effcustomtheta0(i)
      else if (labelfound) then
        ufftheta0 = labeltheta0
      else
        ufftheta0 = 0.0_wp !> no angle target available -- distance restraint only (gated below)
        if (have_metal_r1) write (stdout,'(3x,a,i0,a)') &
        & '**NOTE** atom ',i,': custom r1 given but no theta0 (from override, UFF4MOF label, or '// &
        & '"geometry=") -- distance restraint only, no angle restraint'
      end if
      if (.not.have_metal_r1) write (stdout,'(3x,a,i0,a)') &
      & '**NOTE** atom ',i,': no UFF4MOF r1 (no resolved label, no custom r1) -- donor detection '// &
      & 'falls back to already-bonded atoms only (no long/weak dative-contact window search), '// &
      & 'and no distance restraint is added for any of them, angle restraint (from theta0/"geometry=") only'
      ndonors = 0
      ncarbondonors = 0
      nnodistdonors = 0
      do j = 1,mol%nat
        if (j == i) cycle
        if (isanymetal(mol%at(j))) cycle !> metal-metal contacts not handled here
        is_carbon_donor = .false.
        skip_distance = .false.
!>--- restrict candidate donors to the usual Lewis-basic heteroatoms, OR
!>    a carbon GFN-FF's own topology already bonds to this metal
!>    (carbonyl/cyanide/alkyl/aryl/carbene, incl. oxidative-addition
!>    products) -- a nearby C that is NOT already bonded is essentially
!>    never itself the coordinating atom, so no UFF4MOF-window search
!>    is ever done for carbon the way it is for the other elements
!>    below: there is no tabulated M-C r1 to search a window around,
!>    and real M-C bonds are strong/covalent, so GFN-FF's own topology
!>    (bondmat) is trusted outright instead of guessed at, and no
!>    distance restraint is ever added for carbon regardless of
!>    have_metal_r1.
        select case (i2e(mol%at(j),'nc')) !> i2e(...,'nc') is title-case, e.g. "Cl","Br"
        case ('N','O','S','P','F','Cl','Br','I')
        case ('C')
          if (bondmat(i,j) < 0.5_wp) cycle
          is_carbon_donor = .true.
          skip_distance = .true.
        case default
          cycle
        end select
        if (.not.is_carbon_donor) then
          if (.not.have_metal_r1) then
!>--- no r1 available for this metal at all (see the NOTE above) --
!>    same bondmat-only fallback as carbon, for the same reason: no
!>    bond-length estimate to search a window around, or to restrain
!>    a distance toward, so take whatever GFN-FF's own topology
!>    already considers bonded and restrain its ANGLE only.
            if (bondmat(i,j) < 0.5_wp) cycle
            skip_distance = .true.
          else
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
          end if
        end if
        if (.not.effrestrain(i)) cycle
        if (.not.skip_distance) then
!>--- distance restraint: skipped whenever skip_distance is set --
!>    either a carbon donor (no reliable M-C r1, real covalent bond
!>    doesn't need the help) or any donor found when this metal has no
!>    r1 source at all (see have_metal_r1 above) -- but EITHER WAY
!>    still added to donorlist just below, so the ANGLE-restraint role
!>    assignment (which pair gets 90 vs. 180 deg) correctly accounts
!>    for the site it occupies instead of silently not knowing it's
!>    there (confirmed this omission previously let the matcher force
!>    an UNRELATED pair of real donors into a trans role that didn't
!>    belong to them, because their true trans partner was an
!>    invisible carbon donor).
          call newcons%bondrangeconstraint(i,j,1.25_wp*rij_uff*aatoau,0.0_wp)
          call env%calc%add(newcons)
          nrestrained = nrestrained+1
        else if (is_carbon_donor) then
          ncarbondonors = ncarbondonors+1
        else
          nnodistdonors = nnodistdonors+1
        end if
        ndonors = ndonors+1
        donorlist(ndonors) = j
      end do
      if (ncarbondonors > 0) then
        write (stdout,'(3x,i0,a)') ncarbondonors, &
        & ' carbon-bonded donor(s) included for angle-restraint geometry only (no distance restraint)'
      end if
      if (nnodistdonors > 0) then
        write (stdout,'(3x,i0,a)') nnodistdonors, &
        & ' donor(s) included for angle-restraint geometry only (no r1 for this metal -> no distance restraint)'
      end if
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
!>--- for theta0~90 deg (square-planar CN4, square-pyramidal CN5, or
!>    octahedral CN6 -- the only UFF4MOF geometries with a genuine
!>    antipodal/trans relationship), find the best-fitting perfect
!>    matching (which specific pairs get the 180-deg role) by scoring
!>    EVERY possible matching against the current geometry and keeping
!>    the lowest-error one -- see best_trans_matching below. This
!>    replaces a naive per-pair "closer to 180 or theta0" greedy rule,
!>    which decides each pair independently and can degenerate: for a
!>    perfectly (or nearly) symmetric starting geometry (e.g. exactly
!>    tetrahedral), EVERY pair can independently look "closer to
!>    theta0" at once, greedily assigning ALL pairs target=theta0 and
!>    NONE target=180 -- a geometrically impossible request for a
!>    4-(or 6-)coordinate center (see tetraplatin/'s README entry,
!>    final_tests/). The matching approach instead considers the full
!>    set of pairs jointly, so it always returns a geometrically VALID
!>    4-at-theta0/2-at-180 (or CN6: 12-at-theta0/3-at-180) split, no
!>    matter how symmetric or distorted the starting geometry is. An odd
!>    donor count (e.g. CN=5, square-pyramidal: 4 basal + 1 apical, the
!>    apical site genuinely has no trans partner) is handled by the
!>    same search allowing exactly one donor to be left out of the
!>    "180 deg" role entirely -- see best_trans_matching's own comment.
        if (abs(ufftheta0-90.0_wp) < 1.0_wp) then
          if (allocated(transpartner)) deallocate (transpartner)
          allocate (transpartner(ndonors),source=0)
          call best_trans_matching(donorlist,ndonors,i,ufftheta0,transpartner)
        else if (ndonors == 5.and.abs(ufftheta0-120.0_wp) < 1.0_wp) then
!>--- regime (3) setup: genuine CN=5 trigonal-bipyramidal (e.g. Fe5+3).
!>    Unlike regime (1)'s possibly-multiple trans pairs, TBP has exactly
!>    ONE axial-axial pair among the 5 donors -- find which one best
!>    fits the CURRENT geometry (see best_axial_pair_tbp below).
          call best_axial_pair_tbp(donorlist,ndonors,i,tbp_axial_a,tbp_axial_b)
        end if
        do a = 1,ndonors-1
          do b = a+1,ndonors
            if (abs(ufftheta0-90.0_wp) < 1.0_wp) then
!>--- regime (1): square-planar/square-pyramidal/octahedral -- genuine
!>    cis/trans roles exist
              if (transpartner(a) == b) then
                ang_target = 180.0_wp
              else
                ang_target = ufftheta0
              end if
            else if (ndonors == 5.and.abs(ufftheta0-120.0_wp) < 1.0_wp) then
!>--- regime (3): genuine CN=5 trigonal-bipyramidal -- THREE angle
!>    classes, not two: the chosen axial pair targets 180 (axial-axial),
!>    every pair crossing exactly one axial donor targets 90
!>    (axial-equatorial), and every pair among the remaining three
!>    donors targets theta0=120 (equatorial-equatorial).
              if ((a == tbp_axial_a.and.b == tbp_axial_b).or. &
              &   (a == tbp_axial_b.and.b == tbp_axial_a)) then
                ang_target = 180.0_wp
              else if (a == tbp_axial_a.or.a == tbp_axial_b.or. &
              &        b == tbp_axial_a.or.b == tbp_axial_b) then
                ang_target = 90.0_wp
              else
                ang_target = ufftheta0 !> 120 deg, equatorial-equatorial
              end if
            else
!>--- regime (2): tetrahedral/trigonal/etc -- no trans role exists,
!>    every pair targets theta0 directly
              ang_target = ufftheta0
            end if
            if (effanglefc(i) > 0.0_wp) then
              call newcons%angleconstraint(donorlist(a),i,donorlist(b),ang_target,k=effanglefc(i),tol=angletol)
            else
              call newcons%angleconstraint(donorlist(a),i,donorlist(b),ang_target,tol=angletol)
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

!>--- pi-coordination (eta-bonded ligand) dissociation protection.
!>    Donor detection above never sees these: ordinary sigma-bonded
!>    carbon (CO/CN-/alkyl/aryl/carbene, incl. oxidative-addition
!>    products) IS caught above via bondmat, but an eta-coordinated
!>    ring carbon (Cp, arene) is not, since GFN-FF doesn't register
!>    each individual ring carbon as a single strong bond to the metal
!>    the way it does for a localized sigma donor. Here, for every
!>    [[calculation.nci_metal]] block that gave a pi_atoms list, split
!>    that flat list into separate rings by connectivity (two Cp rings
!>    on the same metal, e.g. nickelocene, are not bonded to each
!>    other), then restrain each resolved metal atom to the LIVE
!>    centroid of each of its own ring's atoms via a centroidrange
!>    constraint (constraints.f90) -- same two-sided logfermi potential
!>    as the normal donor distance restraint above, just measured to a
!>    running average of positions instead of a second fixed atom, so
!>    it moves correctly with the ring through MTD. Distance-only for
!>    now (no angle restraint) -- the point is to stop the WHOLE ring
!>    drifting away, not to pin its internal geometry, which the ring's
!>    own rigidity + GFN-FF's native topology already handle. restrain
!>    defaults to true, same convention as every other field here.
  block
    integer :: npi,a,bb,cc,ncomp,comp_count,sp
    integer,allocatable :: compid(:)
    logical,allocatable :: visited(:)
    integer,allocatable :: stackbuf(:)
    real(wp) :: d0,centroid(3)
    real(wp),parameter :: cc_cutoff = 1.75_wp !> Angstrom, generous C-C/C-X adjacency for ring splitting
    integer :: npirestrained,npirings

    npirestrained = 0
    npirings = 0
    do k = 1,env%nci_metal%nassign
      if (.not.allocated(env%nci_metal%pi_atoms)) cycle
      if (k > size(env%nci_metal%pi_atoms)) cycle
      if (.not.allocated(env%nci_metal%pi_atoms(k)%idx)) cycle
      npi = size(env%nci_metal%pi_atoms(k)%idx)
      if (npi < 1) cycle
      if (.not.env%nci_metal%restrain(k)) cycle

!>--- split this block's flat pi_atoms list into connected rings
      allocate (compid(npi),source=0)
      allocate (visited(npi),source=.false.)
      allocate (stackbuf(npi))
      ncomp = 0
      do a = 1,npi
        if (visited(a)) cycle
        ncomp = ncomp+1
        sp = 1
        stackbuf(1) = a
        visited(a) = .true.
        do while (sp > 0)
          cc = stackbuf(sp)
          sp = sp-1
          compid(cc) = ncomp
          do bb = 1,npi
            if (visited(bb)) cycle
            if (sqrt(sum((mol%xyz(:,env%nci_metal%pi_atoms(k)%idx(cc)) - &
              & mol%xyz(:,env%nci_metal%pi_atoms(k)%idx(bb)))**2))*autoaa < cc_cutoff) then
              visited(bb) = .true.
              sp = sp+1
              stackbuf(sp) = bb
            end if
          end do
        end do
      end do

!>--- one centroidrange constraint per (matched metal atom) x (ring)
      do i = 1,mol%nat
        if (.not.specmatch(i,k)) cycle
        if (.not.isanymetal(mol%at(i))) cycle
        do comp_count = 1,ncomp
          block
            integer,allocatable :: ringatms(:)
            integer :: nr,rr
            nr = count(compid == comp_count)
            allocate (ringatms(nr))
            rr = 0
            do a = 1,npi
              if (compid(a) == comp_count) then
                rr = rr+1
                ringatms(rr) = env%nci_metal%pi_atoms(k)%idx(a)
              end if
            end do
            centroid(:) = 0.0_wp
            do a = 1,nr
              centroid(:) = centroid(:)+mol%xyz(:,ringatms(a))
            end do
            centroid = centroid/real(nr,wp)
            d0 = sqrt(sum((mol%xyz(:,i)-centroid)**2))
            call newcons%centroidrangeconstraint(i,ringatms,1.3_wp*d0,0.6_wp*d0)
            call env%calc%add(newcons)
            npirestrained = npirestrained+1
            npirings = npirings+1
            write (stdout,'(3x,a,i0,a,i0,a,i0,a,f6.2,a)') &
            & 'pi-coordination: atom ',i,' <-> centroid of ',nr,' atoms (ring ',comp_count, &
            & ', d0=',d0*autoaa,' A) restrained'
            deallocate (ringatms)
          end block
        end do
      end do
      deallocate (compid,visited,stackbuf)
    end do
    if (npirings > 0) then
      write (stdout,'(3x,i0,a)') npirings,' pi-coordination centroid restraint(s) added'
    end if
  end block
  write (stdout,*)

!>--- saftey terminations
  call crest_sampling_skip(env,doreturn)
  if (doreturn) return

!>--- sets the MD length according to a flexibility measure. For NCI mode
!>    (which nci_metal is a subtype of) this formula is atom-count-only
!>    (0.10*(rednat+0.1*rednat**2)) -- the sophisticated WBO-based
!>    flexibility measure normal/organic crest runs use is explicitly
!>    skipped for NCI mode, so two metal complexes with the same atom
!>    count but very different numbers of independent donor sites (e.g.
!>    a simple 4-coordinate complex vs. a 6-coordinate one with several
!>    flexible ligands) would otherwise get an IDENTICAL MTD time budget,
!>    regardless of how much more combinatorial arrangement space
!>    (fac/mer-type isomerism, chelate ring conformers, ...) the
!>    higher-coordination one actually has to cover.
  call md_length_setup(env)
!>--- extra nci_metal-specific scaling on top of the above, by total
!>    coordination number summed over all metal centers (a readily
!>    available proxy for "how many independently-arrangeable donor
!>    sites are there to search over" -- CN=4, e.g. a simple square-
!>    planar/tetrahedral center, is the baseline already covered by
!>    everything tested so far and gets NO extra time (scale=1.0);
!>    every CN beyond that adds 15% more MTD time, capped at 3x so a
!>    large multi-metal cluster doesn't blow up unboundedly without a
!>    deliberate user override (length_ps in the toml still wins over
!>    all of this, as it already did before this change).
  ligand_scale = 1.0_wp
  if (total_cn > 4) then
    ligand_scale = min(3.0_wp,1.0_wp+0.15_wp*real(total_cn-4,wp))
  end if
!>--- SEPARATE from donor count: ligand-INTERNAL flexibility (rotatable
!>    bonds, ring systems, ...) -- a low-CN complex can still have a lot
!>    of this (e.g. FeCl3(OEt2)2's ethyl rotation, or a chelate's own
!>    flexible backbone), and the CN-based term above has no way to see
!>    it at all. flexi() is the SAME covalent-flexibility measure normal/
!>    organic crest runs already use (crest/src/flexi.F90) -- reused
!>    here directly rather than reinventing a metric, deliberately
!>    WITHOUT the WBO single point xtbsp() normally pairs it with (that
!>    adds a real extra GFN0 calculation to every run; flexi() still
!>    works without it, just without distinguishing single/double bonds,
!>    which is a fine trade-off to keep this cheap by default).
  allocate (allinclude(mol%nat),source=1)
  call flexi(mol,mol%nat,allinclude,ligflex)
  deallocate (allinclude)
  ligand_scale = ligand_scale*(1.0_wp+0.5_wp*ligflex)
  ligand_scale = min(4.0_wp,ligand_scale)
  if (ligand_scale > 1.01_wp) then
    write (stdout,'(1x,a,i0,a,f5.3,a,f4.2,a)') 'nci_metal: total coordination number ', &
      & total_cn,', ligand flexibility ',ligflex,' -> scaling MTD length by ',ligand_scale,'x'
    !> NOTE: env%tmtd itself is NOT what propagates forward -- md_length_setup
    !> already finalizes env%mdtime (= anint(max(minimum,env%tmtd))*rfac)
    !> internally before returning, and THAT is what env_to_mddat below
    !> reads (via env%mddat%length_ps, unless the user set length_ps
    !> explicitly in the toml, which still takes final precedence over
    !> all of this exactly as before -- see env_to_mddat's own check).
    !> So the scaling has to apply to env%mdtime here, not env%tmtd.
    env%mdtime = env%mdtime*ligand_scale
  end if
!>--- deep_sampling=true: separate, deliberately crude lever from the
!>    CN/flexibility-based scaling above -- directly multiplies the
!>    per-MTD length (3x) and raises the inner mtdloop's iteration cap
!>    (env%Maxrestart, below) regardless of what the flexibility
!>    estimate says, to test whether raw sampling budget alone (no GC)
!>    closes a gap like the one found on Vitamin B12. See classes.f90.
  if (env%deep_sampling) then
    write (stdout,'(1x,a,f6.1,a)') 'deep_sampling=true: multiplying MTD length by 3x (was ', &
    & env%mdtime,' ps)'
    env%mdtime = env%mdtime*3.0_wp
    write (stdout,'(1x,a,i0,a,i0)') 'deep_sampling=true: raising Maxrestart (inner MTD-loop '// &
    & 'iteration cap) from ',env%Maxrestart,' to ',max(env%Maxrestart,15)
    env%Maxrestart = max(env%Maxrestart,15)
  end if
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
!>--- isomer seeding ("seed=true"): decide whether to seed at all, and
!>    enumerate the alternate starting structures ONCE (they don't change
!>    between mainloop/restart iterations, unlike the Kpush/alpha grid
!>    rebuilt each iteration below). Picks the first metal atom with a
!>    resolved seed=true; isomer generation itself is not yet metal-
!>    selective for a heterobimetallic structure (it enumerates ligand
!>    arrangements around ONE given metal_idx) -- the per-atom effseed(:)
!>    is still stored correctly for that future extension, but today's
!>    behavior is "seed around the first such metal, or not at all".
  do_seed = .false.
  seed_metal_idx = 0
  do i = 1,mol%nat
    if (isanymetal(mol%at(i)).and.effseed(i)) then
      do_seed = .true.
      seed_metal_idx = i
      exit
    end if
  end do
  n_seed_structs = 0
  if (do_seed) then
    call generate_isomer_structures(mol%nat,mol%at,mol%xyz,seed_metal_idx, &
    & 100.0_wp,seed_structures,n_seed_structs,seed_geomname,seed_ok)
    if (.not.seed_ok.or.n_seed_structs < 1) then
      write (stdout,'(1x,a)') 'seed=true: no distinct alternate isomer found '// &
      & '(single geometry, or isomer generation not applicable here) -- running unseeded.'
      do_seed = .false.
      n_seed_structs = 0
    else
!>--- break exact symmetry before this candidate is used for anything.
!>    generate_isomer_structures places donors at EXACTLY idealized
!>    angles (90.000000/180.000000 deg, by design -- that's the whole
!>    point of the alignment fix). Confirmed directly (Co(NO2)3(NH3)3's
!>    mer candidate): this exact-180-deg symmetry can trigger a genuine
!>    NaN in GFN-FF's own energy/gradient evaluation (an angle- or
!>    dihedral-dependent term with a removable singularity, e.g.
!>    1/sin(theta), that a real/perturbed structure would never hit,
!>    but an exactly-linear donor arrangement does) -- NOT a convergence
!>    failure, a NaN from the very first energy call, confirmed via
!>    verbose optimizer output. A trivial (<<0.001 Bohr), deterministic,
!>    per-atom offset is enough to clear it without perturbing the
!>    actual chemistry; pre-relax (below) then properly re-optimizes
!>    away this tiny nudge along with everything else.
      block
        integer :: pkk,pia
        real(wp) :: ppert
        do pkk = 1,n_seed_structs
          do pia = 1,mol%nat
            ppert = 1.0e-4_wp*real(mod(pia*37,11)-5,wp)
            seed_structures(1,pia,pkk) = seed_structures(1,pia,pkk)+ppert
            seed_structures(2,pia,pkk) = seed_structures(2,pia,pkk)+ppert*0.7_wp
            seed_structures(3,pia,pkk) = seed_structures(3,pia,pkk)+ppert*1.3_wp
          end do
        end do
      end block
!>--- drop any generated candidate that is a PURE enantiomer of the input
!>    (mirror image only, no genuine donor-arrangement difference) when
!>    distinguish_enantiomers is false (the default): CREGEN's own final
!>    merge will fold it back into the input's entry anyway (iinversion
!>    stays at its normal 0/auto default in that case), so running a full
!>    Kpush/alpha MTD grid on it is pure wasted compute. Genuine
!>    diastereomers (e.g. mer vs fac, or Co(pn)3-type stereocenters) are
!>    NOT mirror images of the input and are never caught by this filter.
      if (.not.env%distinguish_enantiomers) then
        block
          type(coord) :: ref_mol,cand_mol
          real(wp) :: mirror_rmsd
          real(wp),parameter :: enantiomer_rmsd_thr = 0.6_wp !> Bohr (~0.3 A)
          integer :: kk,kk2,io,n_dropped,n_kept,n_invalid
          logical,allocatable :: keep(:)
          real(wp),allocatable :: filtered(:,:,:)
          real(wp),allocatable :: cn_check(:),bond_check(:,:)
          ref_mol = mol
          allocate (keep(n_seed_structs),source=.true.)
          n_dropped = 0
          n_invalid = 0
          do kk = 1,n_seed_structs
            cand_mol = mol
            cand_mol%xyz = seed_structures(:,:,kk)
!>--- self-contained sanity check, entirely on our own side: a
!>    generated candidate with literally no perceivable bond anywhere
!>    is a degenerate/badly-placed geometry (can also break the
!>    comparison machinery below), not a usable seed. Cancel just this
!>    one candidate instead of passing it on.
            if (allocated(cn_check)) deallocate (cn_check)
            if (allocated(bond_check)) deallocate (bond_check)
            call cand_mol%cn_to_bond(cn_check,bond_check,'cov')
            if (.not.any(cn_check > 0.5_wp)) then
              write (stdout,'(3x,a,i0,a)') 'seed candidate ',kk, &
              & ': no perceivable bonds anywhere in the generated geometry -- invalid, cancelling this one'
              keep(kk) = .false.
              n_invalid = n_invalid+1
              cycle
            end if
            mirror_rmsd = irmsd(ref_mol,cand_mol,iinversion=0,io=io)
!>--- irmsd/min_rmsd can return a -1.0 sentinel (its own internal
!>    quaternion-eigensolver failed to converge) WITHOUT setting io
!>    nonzero -- that failure never propagates up as an error. A real
!>    RMSD is never negative, so guard against it explicitly here
!>    instead of trusting it as "a very close match": misreading -1.0
!>    as < enantiomer_rmsd_thr would silently drop every candidate as
!>    a false "enantiomer", exactly the failure mode seen on B12.
            if (io == 0.and.mirror_rmsd >= 0.0_wp.and.mirror_rmsd < enantiomer_rmsd_thr) then
              keep(kk) = .false.
              n_dropped = n_dropped+1
            end if
          end do
          if (n_invalid > 0) then
            write (stdout,'(1x,a,i0,a)') 'cancelled ',n_invalid, &
            & ' invalid (bondless) generated candidate(s) -- see per-candidate notes above.'
          end if
          if (n_dropped+n_invalid > 0) then
            n_kept = n_seed_structs-n_dropped-n_invalid
            if (n_dropped > 0) then
              write (stdout,'(1x,a,i0,a)') 'distinguish_enantiomers=false: dropped ',n_dropped, &
              & ' generated candidate(s) that are pure enantiomers of the input (would merge '// &
              & 'back into it anyway at the final CREGEN step; not worth a full MTD grid).'
            end if
            if (n_kept < 1) then
              do_seed = .false.
              n_seed_structs = 0
              write (stdout,'(1x,a)') 'seed=true: no non-enantiomeric alternate isomer left '// &
              & 'after filtering -- running unseeded.'
            else
              allocate (filtered(3,mol%nat,n_kept))
              kk2 = 0
              do kk = 1,n_seed_structs
                if (keep(kk)) then
                  kk2 = kk2+1
                  filtered(:,:,kk2) = seed_structures(:,:,kk)
                end if
              end do
              call move_alloc(filtered,seed_structures)
              n_seed_structs = n_kept
            end if
          end if
        end block
      end if
      if (do_seed) then
      write (stdout,'(1x,a,i0,a,i0,a)') 'seed=true: found ',n_seed_structs, &
      & ' distinct alternate starting isomer(s) around atom ',seed_metal_idx, &
      & ' ('//trim(seed_geomname)//') -- the full Kpush/alpha grid will run on '// &
      & 'each of them, in addition to the input structure.'
!>--- GFN-FF perceives molecular topology ONCE per structure, from
!>    whichever geometry that structure's own force-field setup first
!>    sees (api_helpers.f90:gfnff_init) -- and that perception can
!>    genuinely differ between seeds that are geometrically distinct
!>    even when they share the exact same molecular graph (confirmed
!>    directly: two independent cis/trans-isomer runs of the same
!>    connectivity landed 75+ kcal/mol apart in absolute energy despite
!>    normal, non-dissociated geometries -- see final_tests/ gaps notes).
!>    Each seed here is deliberately left free to (re-)perceive its own
!>    topology independently -- the physically correct choice for the
!>    MTD *sampling* itself, since forcing one seed's topology onto a
!>    geometrically different seed risks wrong force-field parameters
!>    for that seed's entire exploration, not just a wrong final energy.
!>    The unavoidable cost: GFN-FF energies are NOT guaranteed comparable
!>    across seeds in the resulting ensemble.
      if (any(env%calc%calcs(1:env%calc%ncalculations)%id == jobtype%gfnff)) then
        write (stdout,'(1x,a)') '**NOTE** seed=true with GFN-FF: this ensemble is valid for '// &
        & 'conformer generation (GFN-FF is routinely used this way), but structures from '// &
        & 'different seeds are NOT guaranteed to be on a comparable energy scale -- each '// &
        & 'seed perceives its own topology independently. To get energies you can trust '// &
        & 'for cross-seed comparison, add a second [[calculation.level]] with a '// &
        & 'topology-independent method (e.g. gfn2, the common choice) and refine="sp".'
      end if
      end if
    end if
  end if
  n_total_structs = n_seed_structs+1 !> +1 = the input's own arrangement

!>--- per-seed independent topology: generate_isomer_structures
!>    deliberately EXPANDS donor distances past realistic bond lengths
!>    when repositioning a chelate ligand, to avoid clashes (see
!>    isomer_generate.f90) -- confirmed unsafe to hand straight to the
!>    MTD/reoptimization pipeline as-is (silent donor dissociation
!>    unrestrained; a parallel-reoptimization crash even restrained).
!>
!>    Each block below (1=the original input, 2..n_total_structs=each
!>    generated isomer) gets its OWN independently-initialized
!>    calculator: topology is perceived from the TRUSTED ORIGINAL INPUT
!>    geometry (every block shares the same molecular graph by
!>    construction, so this is a valid topology for all of them), not
!>    from isomer_generate's clash-avoiding expanded geometry -- then
!>    that same, now topology-locked calculator relaxes the block's own
!>    actual starting structure back to something sensible. This is also
!>    the physically correct setup for the MTD *sampling* that follows:
!>    each block explores under its own properly-initialized force
!>    field, rather than all of them being silently forced onto
!>    whichever topology one shared calculator happened to perceive
!>    first (confirmed elsewhere that GFN-FF topology can genuinely
!>    differ between geometrically distinct structures sharing the same
!>    graph -- see the NOTE above). The unavoidable cost remains cross-
!>    block energy comparability -- see that same NOTE and refine="sp".
  if (do_seed) then
    allocate (seedcalcs(n_total_structs))
    block
      type(coord) :: warmmolnew
      real(wp),allocatable :: relaxgrad(:,:)
      real(wp) :: relaxenergy
      integer :: relaxio,kk
      type(constraint),allocatable :: savedcons(:)
      integer :: savednconstraints
      do kk = 1,n_total_structs
        call seedcalcs(kk)%copy(env%calc)
        allocate (relaxgrad(3,mol%nat),source=0.0_wp)
!>--- warm-up: perceive this block's topology from the trusted original
!>    input (same molecular graph for every block)
        call optimize_geometry(mol,warmmolnew,seedcalcs(kk),relaxenergy,relaxgrad, &
        & .false.,.false.,relaxio)
        deallocate (relaxgrad)
        if (kk == 1) cycle !> block 1 = the original input itself, already sensible
!>--- strip this seed's restraint(s) for the pre-relax settling below
!>    ONLY -- confirmed directly to be the actual root cause of the
!>    shallow trap, not a sampling-depth problem: the restraint's
!>    donor-metal-donor angle target is built ONCE, globally, from the
!>    idealized UFF4MOF geometry, and applied identically to every
!>    seed block regardless of isomer type. A generated alternate
!>    isomer's TRUE relaxed geometry can need angles several degrees
!>    off that idealized target (confirmed: up to 6.15 deg on
!>    Co(NO2)3(NH3)3's mer) -- not enough to escape the restraint's
!>    flat-bottom tolerance at the FINAL structure, but enough,
!>    combined across every simultaneously-restrained pair along the
!>    relaxation PATH from the idealized starting candidate, to trap
!>    it in a shallow, restraint-compliant local minimum instead (this
!>    was tested directly: widening the tolerance from 4 to 15 degrees
!>    only partially closed the gap, confirming the issue is path-
!>    dependent resistance, not the final tolerance). Removing the
!>    restraint here only, then restoring it before the main MTD search
!>    below, lets this seed find its own true local minimum first,
!>    exactly like a plain restrain=false run does -- confirmed
!>    directly: this one change alone recovers the true mer minimum
!>    reliably (was: stuck 24 kcal/mol too high every time).
!>    samebinding (below) is the safety net against the one real risk
!>    this removal reintroduces -- losing the donor set entirely during
!>    this brief unrestrained settling -- so it is NOT weakened by this
!>    change: it still runs, after every attempt, before anything is
!>    accepted.
        savednconstraints = seedcalcs(kk)%nconstraints
        if (allocated(seedcalcs(kk)%cons)) then
          call move_alloc(seedcalcs(kk)%cons,savedcons)
        end if
        seedcalcs(kk)%nconstraints = 0
!>--- relax this seed's actual (expanded-donor) geometry using its own
!>    now topology-locked calculator. seed_structures was already
!>    desymmetrized above (tiny offset, just enough to avoid the exact-
!>    symmetry NaN). That alone isn't enough to reliably reach the
!>    TRUE minimum for this donor arrangement, though: a SINGLE
!>    relax attempt starting almost exactly at the idealized candidate
!>    can fall into a shallow, partially-desymmetrized local minimum
!>    instead (confirmed directly on Co(NO2)3(NH3)3's mer candidate --
!>    one attempt converged to a structure keeping only 1 of the 3
!>    trans-donor relationships a genuine mer needs, at an energy even
!>    higher than the OTHER isomer's own basin, instead of the true,
!>    8.6 kcal/mol-lower mer minimum a hand-built mer start finds
!>    directly). Several independent attempts, each with its own
!>    larger, genuinely different random displacement, cheaply
!>    approximate a local basin-hopping step around the idealized
!>    candidate -- keep whichever attempt converges lowest.
        block
          integer,parameter :: nattempts = 20
          type(coord) :: trymol,trymolnew
          type(calcdata) :: trycalc
          real(wp),allocatable :: trygrad(:,:)
          real(wp) :: tryenergy,bestenergy
          integer :: tryio,iatt,ia,ja,idon,ndon,nq,qh
          real(wp) :: dij
          logical :: anyok,samebinding
          logical,allocatable :: origbonded(:),trybonded(:)
          integer,allocatable :: donoratoms(:),queue(:)
          logical,allocatable :: infrag(:,:)
          real(wp) :: axis(3),axislen,angle,ca,sa,rel(3),crossv(3),rotated(3)
!>--- record which atoms are bonded to the metal in the ORIGINAL
!>    (pre-perturbation) candidate, so attempts that are large enough
!>    to change the actual BINDING MODE (not just the rotamer) can be
!>    rejected. Confirmed necessary directly: an earlier, more
!>    aggressive version of a plain-Cartesian-jitter attempt let one
!>    NO2 ligand flip from N-bound (nitro) to O-bound (nitrito) during
!>    relaxation -- a real, lower-energy structure, but a genuinely
!>    different chemical species, not a deeper rotamer of the isomer
!>    this block is supposed to be seeding.
          allocate (origbonded(mol%nat),source=.false.)
          allocate (trybonded(mol%nat),source=.false.)
          do ia = 1,mol%nat
            if (ia == seed_metal_idx) cycle
            dij = norm2(seed_structures(:,ia,kk-1)-seed_structures(:,seed_metal_idx,kk-1))
            if (dij < 3.0_wp*aatoau) origbonded(ia) = .true.
          end do
!>--- identify each donor's own ligand fragment (bondmat reflects the
!>    molecular GRAPH, which is identical for every isomer -- only the
!>    3D arrangement differs -- so it's safe to reuse here even though
!>    it was computed from the original input's own geometry). This
!>    targets the actual missing degree of freedom directly: plain
!>    Cartesian jitter on all atoms either stays in the same shallow
!>    rotamer (too small) or crosses into different chemistry entirely
!>    (too large, confirmed directly) without ever properly sampling
!>    DIFFERENT ROTATIONS of a ligand about its own Co-donor axis --
!>    exactly the freedom a hand-built starting structure gets right
!>    by construction and this generator does not.
          ndon = count(bondmat(1:mol%nat,seed_metal_idx) > 0.5_wp)
          allocate (donoratoms(ndon))
          idon = 0
          do ia = 1,mol%nat
            if (ia == seed_metal_idx) cycle
            if (bondmat(ia,seed_metal_idx) > 0.5_wp) then
              idon = idon+1
              donoratoms(idon) = ia
            end if
          end do
          allocate (infrag(mol%nat,ndon),source=.false.)
          allocate (queue(mol%nat))
          do idon = 1,ndon
            infrag(donoratoms(idon),idon) = .true.
            queue(1) = donoratoms(idon)
            nq = 1
            qh = 1
            do while (qh <= nq)
              ia = queue(qh)
              qh = qh+1
              do ja = 1,mol%nat
                if (ja == seed_metal_idx) cycle !> never cross back through the metal
                if (infrag(ja,idon)) cycle
                if (bondmat(ia,ja) > 0.5_wp) then
                  infrag(ja,idon) = .true.
                  nq = nq+1
                  queue(nq) = ja
                end if
              end do
            end do
          end do
          anyok = .false.
          bestenergy = huge(1.0_wp)
          do iatt = 1,nattempts
            call trycalc%copy(seedcalcs(kk))
            trymol%nat = mol%nat
            trymol%at = mol%at
            trymol%xyz = seed_structures(:,:,kk-1)
            do idon = 1,ndon
              !> deterministic, attempt- and donor-dependent pseudo-
              !> random rotation angle for THIS ligand about its own
              !> Co-donor axis -- the donor atom itself sits ON that
              !> axis (by construction below) and so stays exactly
              !> where it already correctly is; only the rest of the
              !> ligand turns.
              angle = real(mod(idon*331+iatt*977,360),wp)*3.14159265358979_wp/180.0_wp
              axis = trymol%xyz(:,donoratoms(idon))-trymol%xyz(:,seed_metal_idx)
              axislen = norm2(axis)
              if (axislen < 1.0e-6_wp) cycle
              axis = axis/axislen
              ca = cos(angle); sa = sin(angle)
              do ia = 1,mol%nat
                if (.not.infrag(ia,idon)) cycle
                if (ia == donoratoms(idon)) cycle
                rel = trymol%xyz(:,ia)-trymol%xyz(:,donoratoms(idon))
                !> Rodrigues' rotation formula about `axis`
                crossv(1) = axis(2)*rel(3)-axis(3)*rel(2)
                crossv(2) = axis(3)*rel(1)-axis(1)*rel(3)
                crossv(3) = axis(1)*rel(2)-axis(2)*rel(1)
                rotated = rel*ca+crossv*sa+axis*dot_product(axis,rel)*(1.0_wp-ca)
                trymol%xyz(:,ia) = trymol%xyz(:,donoratoms(idon))+rotated
              end do
            end do
            allocate (trygrad(3,mol%nat),source=0.0_wp)
            call optimize_geometry(trymol,trymolnew,trycalc,tryenergy,trygrad, &
            & .false.,.false.,tryio)
            deallocate (trygrad)
            if (tryio == 0) then
              trybonded = .false.
              do ia = 1,mol%nat
                if (ia == seed_metal_idx) cycle
                dij = norm2(trymolnew%xyz(:,ia)-trymolnew%xyz(:,seed_metal_idx))
                if (dij < 3.0_wp*aatoau) trybonded(ia) = .true.
              end do
              samebinding = all(trybonded.eqv.origbonded)
              if (samebinding.and.(.not.anyok.or.tryenergy < bestenergy)) then
                anyok = .true.
                bestenergy = tryenergy
                seed_structures(:,:,kk-1) = trymolnew%xyz
                relaxio = 0
              end if
            end if
          end do
          if (.not.anyok) then
            relaxio = -1
          else
            block
              character(len=64) :: dbgname
              integer :: dbgu,dbgia
              write(dbgname,'(a,i0,a)') 'DEBUG_postbasinhop_seed',kk-1,'.xyz'
              open(newunit=dbgu,file=trim(dbgname),status='replace')
              write(dbgu,'(i0)') mol%nat
              write(dbgu,'(a,f20.10)') 'post-basinhop energy=',bestenergy
              do dbgia=1,mol%nat
                write(dbgu,'(a2,3f20.10)') i2e(mol%at(dbgia),'nc'),seed_structures(:,dbgia,kk-1)*autoaa
              end do
              close(dbgu)
            end block
          end if
        end block
!>--- restore the restraint(s) stripped above, now that this seed's own
!>    true local minimum has been found -- but re-centered on THIS
!>    seed's own actual relaxed angles, not the generic idealized
!>    UFF4MOF target every block started from. Same idea as the
!>    existing per-seed topology warm-up above (each block gets its
!>    own correctly-perceived GFN-FF topology instead of one shared
!>    one) -- here applied to the restraint's angle target instead of
!>    to GFN-FF's topology. This keeps the restraint doing its real
!>    job (stopping the long main MTD search below from drifting into
!>    dissociation or a different donor set) without it fighting the
!>    isomer's own correct geometry the way a fixed idealized target
!>    does. Distance restraints are left untouched -- only the angle
!>    target was shown to cause the problem.
        seedcalcs(kk)%nconstraints = savednconstraints
        if (allocated(savedcons)) then
          block
            integer :: ic,ia1,ia2,ia3
            real(wp) :: v1(3),v2(3),nowang
            do ic = 1,size(savedcons)
              if (savedcons(ic)%type == 2) then !> angle constraint (constraints.f90: angle=2)
                ia1 = savedcons(ic)%atms(1)
                ia2 = savedcons(ic)%atms(2)
                ia3 = savedcons(ic)%atms(3)
                v1 = seed_structures(:,ia1,kk-1)-seed_structures(:,ia2,kk-1)
                v2 = seed_structures(:,ia3,kk-1)-seed_structures(:,ia2,kk-1)
                nowang = acos(max(-1.0_wp,min(1.0_wp, &
                & dot_product(v1,v2)/(norm2(v1)*norm2(v2)))))
                savedcons(ic)%ref(1) = nowang
              end if
            end do
          end block
          call move_alloc(savedcons,seedcalcs(kk)%cons)
        end if
        if (relaxio /= 0) then
          write (stdout,'(1x,a,i0,a)') '**WARNING** seed isomer ',kk-1, &
          & ' failed to pre-relax in any of 5 attempts -- using its unrelaxed '// &
          & '(donor-expanded) geometry as-is; this isomer''s results may be unreliable.'
        end if
      end do
    end block
  end if

!>--- restart detection. Scoped to the unseeded path only: the seeded
!>    path runs a completely different per-block/per-isomer structure
!>    (own trajectories, own calculators, concatenated afterward --
!>    see the do_seed block above) with no equivalent single-stage
!>    checkpoint, and nothing here records which block a kill happened
!>    in. Resuming that correctly would need its own, richer
!>    checkpoint; attempting it with this single-stage one risks
!>    silently mixing blocks instead of just losing the work a restart
!>    is meant to save. do_seed is already fully resolved by this
!>    point (isomer generation + enantiomer filtering above), so this
!>    check is accurate, not a guess based on the toml alone.
  do_restart = .false.
  skip_collect = .false.
  if (do_seed) then
    if (restart_file_exists()) then
      write (stdout,'(1x,a)') '**NOTE** crest.restart found but seed=true '// &
      & 'is active this run -- restart is not yet supported for the seeded '// &
      & 'path, starting fresh (existing checkpoint left untouched).'
    end if
  else if (env%allowrestart.and.restart_file_exists()) then
    call read_restart_log(rdat)
    if (rdat%runtype == crest_nci_metal.and.rdat%stage /= 'done') then
      do_restart = .true.
      call print_restart_info(rdat)
      skip_collect = (rdat%stage == 'post_collect')
    end if
  end if

!>--- Start mainloop
  env%nreset = 0
  start = .true.
! ── apply restart state ───────────────────────────────────────────
  if (do_restart) then
    env%nreset    = rdat%main_iter
    env%elowest   = rdat%elowest
    env%eprivious = rdat%eprivious
    env%nmetadyn  = rdat%nmetadyn
    start = .false.
    call restart_restore_reference(env,rdat%last_file)
    call env%ref%to(mol)
  end if
  MAINLOOP: do
    call printiter
    if (do_restart) then
!>--- restart: preserve .cre_*.xyz files, skip cleanup
      continue
    else if (.not.start) then
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

! ── restart: skip iterations already fully checkpointed ───────────
      if (do_restart) then
        if (rdat%stage == 'mtd_loop'.and.i <= rdat%mtd_iter) cycle mtdloop
        if (rdat%stage == 'mtd_trj'.and.i < rdat%mtd_iter) cycle mtdloop
      end if

      write (stdout,*)
      write (stdout,'(1x,a)') '------------------------------'
      write (stdout,'(1x,a,i0)') 'Meta-Dynamics Iteration ',i
      write (stdout,'(1x,a)') '------------------------------'

      call tim%start(2,'Metadynamics (MTD)')
      if (.not.do_seed.and.do_restart.and.i == rdat%mtd_iter.and. &
      & rdat%stage == 'mtd_trj') then
        write (stdout,'(1x,a,i0,a)') 'Restarting iteration ',i, &
        & ' from existing trajectory/ensemble'
        ensnam = trim(rdat%last_file)
      else if (.not.do_seed) then
        nsim = -1 !>--- enambles automatic MTD setup in init routines
        call crest_search_multimd_init(env,mol,mddat,nsim)
        allocate (mddats(nsim),source=mddat)
        call crest_search_multimd_init2(env,mddats,nsim)
        call crest_search_multimd(env,mol,mddats,nsim)
!>--- a file called crest_dynamics.trj.xyz should have been written
        ensnam = 'crest_dynamics.trj.xyz'
        if (allocated(mddats)) deallocate (mddats)
!>--- checkpoint: trajectory ready, optimization about to start
        call write_restart_log(crest_nci_metal,'mtd_trj',env%nreset,i, &
          &  env%nmetadyn,env%elowest,env%eprivious,ensnam)
      else
        nsim = -1 !>--- enambles automatic MTD setup in init routines
        call crest_search_multimd_init(env,mol,mddat,nsim)
        allocate (mddats(nsim),source=mddat)
        call crest_search_multimd_init2(env,mddats,nsim)
!>--- seeded dispatch: the FULL bias-parameter grid just built above runs
!>    on EVERY starting structure (input + each generated isomer), i.e.
!>    total job count = nsim*n_total_structs, NOT the grid divided across
!>    structures -- confirmed as the intended design (each isomer gets
!>    the same screening the unseeded case would have gotten).
!>
!>    Each block runs as its OWN crest_search_multimd2 call (plural
!>    mols(:), one per replicate -- already used by search_conformers.f90/
!>    search_entropy.f90 for the per-replicate-different-structure
!>    pattern), with env%calc temporarily swapped to that block's own
!>    independently-initialized calculator (seedcalcs(:), see above) --
!>    necessary because crest_search_multimd2 always derives its per-
!>    thread calculators from env%calc directly, with no per-structure
!>    calculator argument of its own. crest_search_multimd2 OVERWRITES
!>    (not appends) crest_dynamics.trj.xyz on every call, so each block's
!>    trajectory is renamed aside immediately and all blocks are
!>    concatenated into the one combined file the rest of the pipeline
!>    expects once every block has run.
        block
          type(calcdata) :: origcalc
          character(len=64) :: blockfile,blockrot
          character(len=:),allocatable :: catcmd
          integer :: catio
          call origcalc%copy(env%calc)
          do iblock = 1,n_total_structs
            block
              type(coord) :: blockmols(nsim)
              type(mddata) :: blockmddats(nsim)
              integer :: jj
              do jj = 1,nsim
                blockmols(jj) = mol !> copies at/nat/lat/etc.; xyz overwritten next
                if (iblock == 1) then
                  blockmols(jj)%xyz = mol%xyz
                else
                  blockmols(jj)%xyz = seed_structures(:,:,iblock-1)
                end if
                blockmddats(jj) = mddats(jj)
              end do
              call env%calc%copy(seedcalcs(iblock))
!>--- NOTE: seedcalcs(iblock)'s angle restraint (if iblock>1) was
!>    already re-centered on this seed's own true relaxed angles above
!>    (see the pre-relax restore block), not left at the generic
!>    idealized UFF4MOF target -- so it is kept active here, same as
!>    any other block, rather than stripped. An earlier version of
!>    this fix stripped it entirely for the whole MTD+reopt stage
!>    (confirmed necessary at the time: restoring the UNMODIFIED,
!>    idealized-target restraint here pulled the seed straight back to
!>    its old shallow trap, -3.2574 Ha, matching the pre-fix value to
!>    4 decimals, even though pre-relax had already found the true
!>    minimum) -- re-centering the target instead keeps the restraint's
!>    real job (stopping this block's MTD from drifting into
!>    dissociation or a different donor set) without that fight.
              call crest_search_multimd2(env,blockmols,blockmddats,nsim)
            end block
            write (blockfile,'(a,i0,a)') 'crest_dynamics_seedblock',iblock,'.trj.xyz'
            call rename('crest_dynamics.trj.xyz',trim(blockfile))
!>--- reoptimize THIS block's own raw MTD trajectory NOW, while
!>    env%calc is STILL this block's own seed-specific calculator (set
!>    just above) -- not after it gets restored to the shared/original
!>    one below. This is the fix for a real bug found via a direct
!>    debug dump of this exact pipeline: a seed's genuinely distinct
!>    local minimum (a mer-type donor arrangement, correctly generated
!>    by isomer_generate.f90 and correctly surviving the per-seed
!>    pre-relax above: 1/3 raw-generated candidates stayed mer) still
!>    silently vanished (0/5360 in the final ensemble) under the OLD
!>    code, which deferred ALL reoptimization to a single generic pass
!>    running AFTER env%calc was already restored to the ORIGINAL
!>    (differently-perceived) topology -- re-relaxing every frame,
!>    regardless of origin block, under force-field parameters fixed
!>    from a different arrangement that had no stable minimum there.
!>    CREGEN's own sort/dedup pass never touches env%calc (no geometry
!>    optimization happens inside it, confirmed in cregen.f90/
!>    confscript2_misc.f90), so merging already-correctly-reoptimized
!>    per-block results below (the new order) is safe where merging the
!>    RAW trajectories before any reoptimization (the old order) was not.
            call optlev_to_multilev(env%optlev,multilevel)
            call crest_multilevel_oloop(env,trim(blockfile),multilevel,0,i)
            if (env%iostatus_meta .ne. 0) return
!>--- ALSO apply the final very-tight optimization level HERE, per
!>    block, while env%calc is STILL this block's own correctly-
!>    perceived calculator (swapped in above). The separate final-
!>    ensemble optimization further down (crest_multilevel_wrap, run
!>    once after all blocks are merged) is SKIPPED for the seeded path
!>    for exactly this reason: by the time that call would run, env%calc
!>    has been restored to a single shared calculator (origcalc, below),
!>    and re-optimizing a mixed-origin pool under one topology would
!>    silently reintroduce the same per-seed-topology bug this whole
!>    per-block mechanism exists to avoid -- just one stage later than
!>    the original fix caught it. Confirmed directly: a seed's mer-type
!>    donor arrangement, stable and correctly reoptimized at this point,
!>    was still lost (0% survival in the final ensemble) before this
!>    change, despite surviving in the pooled pre-final-optimization set
!>    (43.3% mer) -- i.e. lost specifically during that later, shared-
!>    calculator pass.
            call checkname_xyz(crefile,atmp,str)
            call crest_multilevel_wrap(env,trim(atmp),0)
            if (env%iostatus_meta .ne. 0) return
            call checkname_xyz(crefile,atmp,str)
            write (blockrot,'(a,i0,a)') 'crest_rotamers_seedblock',iblock,'.xyz'
            call rename(trim(atmp),trim(blockrot))
          end do
          call env%calc%copy(origcalc)
!>--- merge all blocks' OWN already-correctly-reoptimized conformer
!>    sets (not the raw trajectories -- reoptimization already happened
!>    per block, above) into the next free crefile slot, where the
!>    unchanged code below expects to find this iteration's result.
!>    Via a shell 'cat', not a manual line-by-line Fortran copy: a
!>    fixed-length read buffer silently truncates any source line
!>    longer than its declared length with no I/O error raised, which
!>    corrupted larger (many-atom) ensembles' extxyz frames outright
!>    (confirmed directly: a 79-atom system's recombined trajectory
!>    produced "unexpected line format" parse errors downstream) --
!>    'cat' has no such limit.
          catcmd = 'cat'
          do iblock = 1,n_total_structs
            write (blockrot,'(a,i0,a)') 'crest_rotamers_seedblock',iblock,'.xyz'
            catcmd = catcmd//' '//trim(blockrot)
          end do
          call checkname_xyz(crefile,atmp,str)
          catcmd = catcmd//' > '//trim(str)
          call execute_command_line(catcmd,exitstat=catio)
!>--- keep each block's OWN already-correctly-reoptimized rotamers
!>    file as a separate, permanent output file, instead of deleting
!>    it once merged into the pooled one above. Requested directly:
!>    the pooled file's single shared energy window can make an
!>    entire seeded isomer type vanish (correctly, by its own
!>    definition -- see gaps.html) if it happens to sit further than
!>    ewin above whatever the pool's lowest structure turns out to be,
!>    which conflates "is this a reasonable structure of THIS isomer"
!>    with "did this isomer survive at all" -- two different
!>    questions. These per-seed files answer the first question alone.
!>    NOTE: these are ROTAMERS files (this block's own optimized pool,
!>    pre-dedup), the same meaning "rotamers" has everywhere else in
!>    CREST -- NOT yet put through CREGEN's RMSD-based duplicate
!>    merge+window the way the final crest_conformers.xyz is (that
!>    dedup only happens once, on the pooled file, below). A real
!>    per-seed CREGEN pass (so each file reports genuinely unique
!>    conformers of just that isomer, windowed against its own
!>    minimum) would be more useful still, but needs its own care to
!>    avoid disturbing the pooled merge this same cregen_interface
!>    machinery also has to do right after -- left for a follow-up,
!>    not done here. Block 1 is the original input's own isomer;
!>    blocks 2+ are the generated alternates, numbered to match
!>    seed_structures(:,:,iblock-1)'s own indexing elsewhere in this
!>    file.
!>    NOTE on the name chosen below: deliberately NOT "crest_rotamers_*"
!>    -- V2terminating (cleanup.f90, shared code, runs at the very end
!>    of every crest invocation) does an unconditional wildcard
!>    rmrfw('crest_rotamers_'), which would silently delete these again
!>    right after writing them if they matched that prefix (confirmed
!>    directly: an earlier version of this using that prefix produced
!>    the log message and the files, then they were gone by the time
!>    the run actually finished). "nci_metal_seed*" matches none of
!>    V2terminating's or V2cleanup's own wildcard patterns.
          do iblock = 1,n_total_structs
            write (blockrot,'(a,i0,a)') 'crest_rotamers_seedblock',iblock,'.xyz'
            write (blockfile,'(a,i0,a)') 'nci_metal_seed',iblock-1,'_rotamers.xyz'
            call rename(trim(blockrot),trim(blockfile))
          end do
          write (stdout,'(/,1x,a)') 'Per-seed rotamer pools kept as separate files (each block''s '// &
          & 'own optimized structures, NOT yet CREGEN-deduplicated, so an isomer can never be '// &
          & 'silently absent from every output at once):'
          do iblock = 1,n_total_structs
            write (blockfile,'(a,i0,a)') 'nci_metal_seed',iblock-1,'_rotamers.xyz'
            write (stdout,'(3x,a,i0,a,a)') 'seed ',iblock-1,': ',trim(blockfile)
          end do
        end block
      end if
      call tim%stop(2)
!>--- deallocate for next iteration
      if (allocated(mddats)) deallocate (mddats)

!==========================================================!
!>--- Reoptimization of trajectories (the seeded case already did
!>    this, per block, under each block's own correct calculator,
!>    above -- doing it again here with one shared calculator would
!>    silently re-introduce the exact bug that per-block step fixes)
      call tim%start(3,'Geometry optimization')
      if (.not.do_seed) then
!>--- a file called crest_dynamics.trj.xyz should have been written
        ensnam = 'crest_dynamics.trj.xyz'
        call optlev_to_multilev(env%optlev,multilevel)
        call crest_multilevel_oloop(env,ensnam,multilevel,i,i)
        if (env%iostatus_meta .ne. 0) return
      end if
      call tim%stop(3)

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
!>--- checkpoint after this MTD iteration (nmetadyn already updated above)
        if (.not.do_seed) call write_restart_log(crest_nci_metal,'mtd_loop', &
          &  env%nreset,i,env%nmetadyn,env%elowest,env%eprivious,trim(str))
!>-- and always do two cycles of MTDs
        cycle mtdloop
      end if
!>--- checkpoint after this MTD iteration
      if (.not.do_seed) call write_restart_log(crest_nci_metal,'mtd_loop', &
        &  env%nreset,i,env%nmetadyn,env%elowest,env%eprivious,trim(str))
!=========================================================!
!>--- Check for lowest energy
      call elowcheck(lower,env)
      if (.not.lower) then
        exit mtdloop
      end if
    end do mtdloop
    do_restart = .false.
!=========================================================!
!>--- collect all ensembles from mtdloop and merge
    write (stdout,*)
    write (stdout,'(''========================================'')')
    write (stdout,'(''           MTD Simulations done         '')')
    write (stdout,'(''========================================'')')
    write (stdout,'(1x,''Collecting ensmbles.'')')
    if (skip_collect) then
!>--- post_collect restart: collectcre/newcregen already ran, reuse last file
      inquire (file=trim(rdat%last_file),exist=fex)
      if (.not.fex) then
        write (stdout,'(/,a)') '**ERROR** restart ensemble not found: '// &
          & trim(rdat%last_file)
        write (stdout,'(a,/)') ' Delete crest.restart and rerun from scratch.'
        call creststop(status_safety)
      end if
      atmp = trim(rdat%last_file)
      write (stdout,'(1x,a,a)') 'Restarting from ensemble: ',trim(atmp)
      skip_collect = .false.
    else
!>-- collecting all ensembles saved as ".cre_*.xyz"
      call collectcre(env)
!>DEBUG_DUMP: the TRUE combined input (both mtdloop iterations' .cre_N
!>    backups appended together by collectcre) going into the one real
!>    final newcregen call -- compare this directly against the final
!>    output dump below to see exactly what that call does.
      block
        character(len=256) :: dbgcmd
        character(len=80) :: dbg_atmp,dbg_btmp
        call checkname_xyz(crefile,dbg_atmp,dbg_btmp)
        dbgcmd = 'cp '//trim(dbg_atmp)//' DEBUG_FINALMERGE_INPUT.xyz'
        call execute_command_line(trim(dbgcmd))
      end block
      call newcregen(env,0)
      call checkname_xyz(crefile,atmp,btmp)
!>DEBUG_DUMP: the TRUE final output, right after the one real newcregen call
      block
        character(len=256) :: dbgcmd
        dbgcmd = 'cp '//trim(atmp)//' DEBUG_FINALMERGE_OUTPUT.xyz'
        call execute_command_line(trim(dbgcmd))
      end block
!>--- checkpoint after collection and CREGEN (unseeded path only, same
!>    scoping reason as the detection block above)
      if (.not.do_seed) call write_restart_log(crest_nci_metal,'post_collect', &
        &  env%nreset,0,env%nmetadyn,env%elowest,env%eprivious,trim(atmp))
    end if
!>--- remaining number of structures
    call remaining_in(atmp,env%ewin,nallout)

!=========================================================!
!>--- (optional) Perform GC step -- opt-in via full_search=true (see
!>    classes.f90 for why this defaults off). Same mechanism as iMTD-
!>    GC's own GC block (crest_search_imtdgc, search_conformers.f90):
!>    crest_newcross3 crosses+optimizes+appends new structures into the
!>    running crefile, confg_chk3 re-merges via CREGEN. Unlike iMTD-GC,
!>    a new lowest-energy structure found here does NOT restart the
!>    outer search around it -- nci_metal's MAINLOOP is single-pass by
!>    design (no equivalent "cycle MAINLOOP" path exists here) -- the
!>    crossed structures are simply folded into this same final
!>    ensemble, which remaining_in below re-counts.
    if (env%full_search.and.env%performCross) then
      write (stdout,*)
      write (stdout,'(1x,a)') 'full_search=true: running genetic-crossing (GC) step...'
      call tim%start(5,'Genetic crossing (GC)')
      call crest_newcross3(env)
      call tim%stop(5)
      if (env%iostatus_meta /= 0) return
      call confg_chk3(env)
      call elowcheck(lower,env)
      if (lower) write (stdout,'(1x,a)') 'GC found a new lowest-energy structure.'
      call checkname_xyz(crefile,atmp,btmp)
      call remaining_in(atmp,env%ewin,nallout)
    end if

!>--- repeat the seed+GFN-FF cross-seed-energy warning here, right next
!>    to the final ensemble this run just produced -- the first copy
!>    (printed near the top, before any MTD output) is easy to scroll
!>    past in a long run; this one lands right where the user is
!>    actually looking, at the final conformer count.
    if (do_seed) then
      if (any(env%calc%calcs(1:env%calc%ncalculations)%id == jobtype%gfnff)) then
        write (stdout,'(1x,a,i0,a)') '**WARNING** the ',nallout, &
        & ' conformer(s) above came from seed=true with GFN-FF: energies are '// &
        & 'NOT guaranteed comparable across seeds (each perceived its own '// &
        & 'topology independently). Before ranking or comparing these '// &
        & 'energies across seeds, add a second [[calculation.level]] with a '// &
        & 'topology-independent method and refine="sp".'
      end if
    end if

!==========================================================!
!>--- exit mainloop
    exit MAINLOOP
  end do MAINLOOP

!==========================================================!
!>--- final ensemble optimization: SKIPPED for the seeded path, where
!>    every block already went through this same very-tight level
!>    individually, under its own correctly-perceived calculator (see
!>    the per-block call above). Running it again here, on the merged,
!>    mixed-origin pool, would force every structure through whichever
!>    ONE calculator happens to be active (env%calc was restored to
!>    origcalc -- the original/block-1 topology -- once the per-block
!>    loop finished), silently reintroducing the per-seed-topology bug
!>    for every OTHER block's structures.
  if (.not.do_seed) then
    write (stdout,'(/)')
    write (stdout,'(3x,''================================================'')')
    write (stdout,'(3x,''|           Final Geometry Optimization        |'')')
    write (stdout,'(3x,''================================================'')')
    call tim%start(3,'Geometry optimization')
    call checkname_xyz(crefile,atmp,str)
    call crest_multilevel_wrap(env,trim(atmp),0)
    call tim%stop(3)
    if (env%iostatus_meta .ne. 0) return
  end if

!==========================================================!
!>--- run is complete: drop the restart checkpoint (harmless no-op if
!>    the seeded path never wrote one, or if restart was never used)
  call delete_restart_log()

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
!> Finds the best-fitting perfect matching of DONORS(1:NDON) (around
!> metal atom CENTER) into NDON/2 disjoint "trans-axis" pairs, i.e.
!> which SPECIFIC pairs should be restrained toward 180 deg while every
!> other pair targets THETA0 -- for regime (1) (theta0~90 deg,
!> square-planar/octahedral) restraints, replacing a naive per-pair
!> greedy rule that can degenerate for a symmetric starting geometry
!> (see the comment above this routine's call site, and tetraplatin/'s
!> README entry in final_tests/). Exhaustively scores EVERY perfect
!> matching (there are (NDON-1)!! of them -- 3 for NDON=4, 15 for
!> NDON=6, still trivial to enumerate once per metal atom at restraint
!> setup) against the CURRENT geometry and keeps the lowest-error one.
!> BESTPARTNER(k)=m (and BESTPARTNER(m)=k) means donors k and m (LOCAL
!> indices into DONORS(:), i.e. 1..NDON, not global atom indices) are
!> the winning matching's trans-axis partners. Only called for even
!> NDON -- odd counts have no perfect matching and fall back to the old
!> greedy rule at the call site.
!========================================================================================!
!>--- for an ODD donor count (e.g. CN=5, square-pyramidal: 4 basal +
!>    1 apical), one donor genuinely has no trans partner at all (the
!>    apical site is cis, ~90 deg, to all 4 basal ones in the idealized
!>    template -- see get_geometry_sites('square_pyramidal') in
!>    isomer_enumerate.f90). A perfect matching (every donor paired)
!>    does not exist for an odd count, but a matching that pairs ALL
!>    BUT ONE still does, and is exactly what a real 5-coordinate
!>    square-pyramidal center needs: 2 trans (180 deg) pairs among 4 of
!>    the 5 donors, the 5th left out of the "180" role entirely (every
!>    pair involving it falls to theta0 by the normal else-branch in
!>    matching_score below, with no special-casing needed there).
!>    nleftover is 0 for even ndon (the original exact behavior,
!>    unchanged) or 1 for odd ndon; enumerate_matchings tries leaving
!>    out each donor in turn (not just one fixed guess) and keeps
!>    whichever leftover + pairing combination fits the real geometry
!>    best, same spirit as the even-count case already did for which
!>    SPECIFIC pairs get the trans role.
  subroutine best_trans_matching(donors,ndon,center,theta0,bestpartner)
    implicit none
    integer,intent(in) :: donors(:),ndon,center
    real(wp),intent(in) :: theta0
    integer,intent(out) :: bestpartner(:)
    integer :: partner(ndon)
    real(wp) :: bestscore
    integer :: nleftover
    partner = 0
    bestpartner(1:ndon) = 0
    bestscore = huge(1.0_wp)
    nleftover = mod(ndon,2)
    call enumerate_matchings(donors,ndon,center,theta0,partner,bestpartner,bestscore,nleftover)
  end subroutine best_trans_matching

  recursive subroutine enumerate_matchings(donors,ndon,center,theta0,partner,bestpartner,bestscore,nleftover)
    implicit none
    integer,intent(in) :: donors(:),ndon,center
    real(wp),intent(in) :: theta0
    integer,intent(inout) :: partner(:)
    integer,intent(inout) :: bestpartner(:)
    real(wp),intent(inout) :: bestscore
    integer,intent(in) :: nleftover !> 0 or 1: how many donors may still go unpaired
    integer :: k,m
    real(wp) :: score
!>--- find the first still-undecided (0) local donor index
    k = 1
    do while (k <= ndon)
      if (partner(k) == 0) exit
      k = k+1
    end do
    if (k > ndon) then
!>--- every donor paired or deliberately left over -> complete
!>    matching, score it (a leftover donor's partner stays 0, which
!>    matching_score's "partner(aa)==bb" test never satisfies, so its
!>    pairs fall to theta0 automatically, same as any non-trans pair)
      score = matching_score(donors,ndon,center,theta0,partner)
      if (score < bestscore) then
        bestscore = score
        bestpartner(1:ndon) = partner(1:ndon)
      end if
      return
    end if
!>--- try pairing donor k with every still-undecided donor after it
    do m = k+1,ndon
      if (partner(m) /= 0) cycle
      partner(k) = m
      partner(m) = k
      call enumerate_matchings(donors,ndon,center,theta0,partner,bestpartner,bestscore,nleftover)
      partner(k) = 0
      partner(m) = 0
    end do
!>--- odd ndon only, and only if this leftover slot hasn't been used
!>    yet elsewhere in this branch: also try leaving k itself unpaired
!>    (sentinel -1, so it's never revisited by the "first undecided"
!>    scan above, which only looks for exactly 0)
    if (nleftover > 0) then
      partner(k) = -1
      call enumerate_matchings(donors,ndon,center,theta0,partner,bestpartner,bestscore,nleftover-1)
      partner(k) = 0
    end if
  end subroutine enumerate_matchings

!> Sum of squared deviations between ACTUAL donor-metal-donor angles and
!> the targets implied by one candidate matching (180 deg for pairs the
!> matching pairs up, THETA0 for every other pair) -- lower is a better
!> fit to that specific matching. Used by enumerate_matchings above.
  real(wp) function matching_score(donors,ndon,center,theta0,partner)
    implicit none
    integer,intent(in) :: donors(:),ndon,center,partner(:)
    real(wp),intent(in) :: theta0
    integer :: aa,bb
    real(wp) :: vv1(3),vv2(3),angnow,targetang,dev
    matching_score = 0.0_wp
    do aa = 1,ndon-1
      do bb = aa+1,ndon
        vv1 = mol%xyz(:,donors(aa))-mol%xyz(:,center)
        vv2 = mol%xyz(:,donors(bb))-mol%xyz(:,center)
        angnow = acos(max(-1.0_wp,min(1.0_wp,dot_product(vv1,vv2)/(norm2(vv1)*norm2(vv2)))))*(180.0_wp/pi)
        if (partner(aa) == bb) then
          targetang = 180.0_wp
        else
          targetang = theta0
        end if
        dev = angnow-targetang
        matching_score = matching_score+dev*dev
      end do
    end do
  end function matching_score

!========================================================================================!
!> Finds, for a genuine CN=5 trigonal-bipyramidal center (regime (3)),
!> which SPECIFIC pair of DONORS(1:NDON) (NDON==5) should hold the
!> axial-axial (180 deg) role, by exhaustively scoring every one of the
!> C(5,2)=10 possible pairs against the CURRENT geometry (the other 3
!> donors are then equatorial: eq-eq pairs target 120 deg, ax-eq cross
!> pairs target 90 deg -- see tbp_axial_score below) and keeping the
!> lowest-error choice. Unlike regime (1)'s best_trans_matching, there
!> is only ONE pair to choose here (not a full matching over all
!> donors), so a plain double loop over the 10 candidate pairs
!> suffices -- no recursive enumeration needed. BESTA/BESTB are LOCAL
!> indices into DONORS(:), i.e. 1..NDON.
!========================================================================================!
  subroutine best_axial_pair_tbp(donors,ndon,center,besta,bestb)
    implicit none
    integer,intent(in) :: donors(:),ndon,center
    integer,intent(out) :: besta,bestb
    integer :: aa,bb
    real(wp) :: score,bestscore
    besta = 1
    bestb = 2
    bestscore = huge(1.0_wp)
    do aa = 1,ndon-1
      do bb = aa+1,ndon
        score = tbp_axial_score(donors,ndon,center,aa,bb)
        if (score < bestscore) then
          bestscore = score
          besta = aa
          bestb = bb
        end if
      end do
    end do
  end subroutine best_axial_pair_tbp

!> Sum of squared deviations between ACTUAL donor-metal-donor angles
!> and the three-class TBP targets (180 for the candidate axial pair
!> AXA/AXB, 90 for every pair involving exactly one of them, 120 for
!> every pair among the remaining three "equatorial" donors) implied
!> by one candidate axial-pair choice. Used by best_axial_pair_tbp.
  real(wp) function tbp_axial_score(donors,ndon,center,axa,axb)
    implicit none
    integer,intent(in) :: donors(:),ndon,center,axa,axb
    integer :: aa,bb
    real(wp) :: vv1(3),vv2(3),angnow,targetang,dev
    tbp_axial_score = 0.0_wp
    do aa = 1,ndon-1
      do bb = aa+1,ndon
        vv1 = mol%xyz(:,donors(aa))-mol%xyz(:,center)
        vv2 = mol%xyz(:,donors(bb))-mol%xyz(:,center)
        angnow = acos(max(-1.0_wp,min(1.0_wp,dot_product(vv1,vv2)/(norm2(vv1)*norm2(vv2)))))*(180.0_wp/pi)
        if ((aa == axa.and.bb == axb).or.(aa == axb.and.bb == axa)) then
          targetang = 180.0_wp
        else if (aa == axa.or.aa == axb.or.bb == axa.or.bb == axb) then
          targetang = 90.0_wp
        else
          targetang = 120.0_wp
        end if
        dev = angnow-targetang
        tbp_axial_score = tbp_axial_score+dev*dev
      end do
    end do
  end function tbp_axial_score
!========================================================================================!
!> Prints, on the current output line (advance='no' style, matching the
!> per-atom "UFF4MOF assignment" printout this is called from), the list
!> of UFF4MOF labels actually tabulated for ELEM -- shown right next to
!> a rejected cross-element/typo'd uff4mof= override so the user sees
!> immediately what they probably meant to type, without having to
!> scroll back up to the earlier "Available UFF4MOF atom types" block.
!========================================================================================!
  subroutine print_uff_typo_hint(elem,opts,nopts)
    implicit none
    character(len=*),intent(in) :: elem
    character(len=8),intent(in) :: opts(:)
    integer,intent(in) :: nopts
    integer :: kk
    if (nopts > 0) then
      write (stdout,'(a,a,a)',advance='no') '. Types tabulated for ',trim(elem),': '
      do kk = 1,nopts
        if (kk > 1) write (stdout,'(a)',advance='no') '/'
        write (stdout,'(a)',advance='no') trim(opts(kk))
      end do
    else
      write (stdout,'(a,a,a)',advance='no') '. No UFF4MOF types are tabulated for ',trim(elem),' at all'
    end if
  end subroutine print_uff_typo_hint
!========================================================================================!
end subroutine crest_search_newnci_metal
!========================================================================================!
