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

!========================================================================================!
!========================================================================================!
module nci_metal_settings_mod
!************************************************************
!* Settings container for the "nci_metal" runtype.
!*
!* Holds the user's choice of UFF4MOF atom type (see uff4mof.f90
!* for the parameter table), set via one or more
!* [[calculation.nci_metal]] TOML blocks. "atom" may be a single
!* 1-based index, a list/range (anything strucrd's get_atlist
!* understands, e.g. "1-4"), or an element symbol (e.g. "Zn",
!* meaning "every Zn atom"), e.g.:
!*
!*   [[calculation.nci_metal]]
!*   atom = 1
!*   uff4mof = "Fe6+2"
!*
!*   [[calculation.nci_metal]]
!*   atom = "Zn"
!*   uff4mof = "Zn6+2"
!*
!*   [[calculation.nci_metal]]
!*   atom = "335-338"
!*   uff4mof = "Fe3+2"
!*
!* The actual structure is not necessarily loaded yet at the point
!* [[calculation.nci_metal]] blocks are parsed (see parse_calcdata.f90),
!* so specs are only stored here as text and resolved to concrete
!* atom IDs later, once the structure is available (see
!* crest_search_newnci_metal). Where multiple specs match the same
!* atom, the one added LAST wins (later blocks override earlier ones).
!*
!* This is a deliberately *minimal* module (only depends on
!* crest_parameters for the working-precision kind wp) so an instance
!* can be stored directly on CREST's systemdata object, analogous to
!* ttconf_settings_mod.
!************************************************************
  use crest_parameters,only:wp
  implicit none
  private

  public :: nci_metal_settings

  !> one block's pi_atoms list (ragged storage: each block may give a
  !> different number of ring/diene carbons, or none at all)
  type :: pi_atomlist
    integer,allocatable :: idx(:)
  end type pi_atomlist

  type :: nci_metal_settings
    integer :: nassign = 0
    character(len=64),allocatable :: spec(:) !> raw atom spec, e.g. "1", "1-4", "Zn"
    character(len=8),allocatable :: uff(:)   !> user-chosen UFF4MOF label
    !> whether the resolved coordination sphere for this spec should be
    !> held together by a distance restraint during MTD/optimization
    !> (see nci_metal_bonds in search_newnci_metal.f90); defaults to
    !> .true. since the point of typing an atom at all is usually to
    !> stop its (possibly weak/long) donor bonds from drifting apart
    !> under the MTD bias potential -- set to .false. per-spec to type
    !> an atom for reporting/auto-detection only, without restraining it.
    logical,allocatable :: restrain(:)
    !> if true (and no explicit uff(:) label was given -- see add_spec),
    !> resolve a CN-ambiguous element (e.g. Cu4+2 vs Cu3f2, both CN4) by
    !> comparing actual donor angles against each candidate's theta0,
    !> instead of leaving it unresolved. See search_newnci_metal.f90.
    logical,allocatable :: auto(:)
    !> if true, an explicit uff(:) label is accepted even if it belongs
    !> to a DIFFERENT element's UFF4MOF namespace (e.g. uff4mof="Pt4+2"
    !> on a Ni atom) -- deliberately borrowing another metal's tabulated
    !> bond radius (r1) and ideal angle (theta0) for creative/exploratory
    !> use. Without this, such a mismatched label is rejected with a
    !> warning and falls back to auto-detection (see search_newnci_metal.f90).
    logical,allocatable :: force(:)
    !> user-settable angle-restraint force constant (Hartree/rad^2),
    !> passed through to angleconstraint's optional k -- default (when
    !> <=0, i.e. not set) is constraints.f90's own fcdefault=0.01, which
    !> is strong enough to HOLD a geometry already close to its target
    !> but far too weak to actively RELOCATE a structure to a
    !> substantially different target angle (e.g. forcing tetrahedral
    !> from an already-square-planar start) -- confirmed empirically:
    !> a 30x longer MTD did not help. Set explicitly, much higher than
    !> 0.01, to actually force a large geometry change; see
    !> search_newnci_metal.f90.
    real(wp),allocatable :: angle_fc(:)
    !> fully user-supplied metal-ligand bond radius (r1, Angstrom) and/or
    !> ideal donor-metal-donor angle (theta0, degrees), bypassing the
    !> UFF4MOF table lookup for that specific quantity entirely. Each is
    !> independent: setting only r1 keeps the auto-detected/tabulated
    !> theta0 (and vice versa), setting both fully replaces the table for
    !> this atom -- useful for a metal/oxidation-state/geometry that
    !> isn't tabulated in UFF4MOF at all, or a deliberately non-standard
    !> (distorted/asymmetric real) target that doesn't match any entry.
    !> Sentinel <=0 (not set) means "use whatever the normal UFF4MOF
    !> resolution -- label/auto/force -- would have produced." Unlike
    !> force=true (which still borrows a real, internally-consistent
    !> tabulated entry), these values are NOT chemically validated by
    !> CREST at all -- the user is fully responsible for physical
    !> sanity. See search_newnci_metal.f90.
    real(wp),allocatable :: r1(:)
    real(wp),allocatable :: theta0(:)
    !> if true, enumerate every symmetry-distinct way this metal's own
    !> ligands could be arranged (isomer_generate_mod::generate_isomer_
    !> structures) and run the FULL existing Kpush/alpha bias-potential
    !> grid on every one of them (not a divided/split grid -- total job
    !> count = grid_size * n_isomers), in addition to the input's own
    !> arrangement, via crest_search_multimd2 instead of the single-
    !> structure crest_search_multimd. Defaults to .true.: the point of
    !> typing a metal at all is usually also to not miss an isomer a
    !> single unconstrained start never crosses to (see e.g. ZnComplex-
    !> seeded/ in final_tests, a fac start that never finds mer on its
    !> own). Set to .false. per-spec for a single plain run, no isomer
    !> enumeration -- see search_newnci_metal.f90.
    !>
    !> GFN-FF caveat: each seed (the input's own arrangement and every
    !> generated isomer) gets its OWN independently (re-)initialized
    !> force-field topology, warmed up from the trusted original input
    !> geometry rather than isomer_generate's clash-avoiding expanded
    !> one -- the physically correct setup for the MTD sampling itself
    !> (confirmed directly: GFN-FF's one-time-per-structure topology
    !> perception can genuinely differ between geometrically distinct
    !> seeds sharing the same molecular graph, e.g. two separate cis/
    !> trans-isomer runs landed 75+ kcal/mol apart despite normal,
    !> non-dissociated geometries). The unavoidable cost: raw GFN-FF
    !> energies are NOT guaranteed comparable across seeds in the
    !> resulting ensemble -- valid for conformer generation (GFN-FF is
    !> routinely used this way) but not for cross-seed energy ranking.
    !> For energies you can trust, add a second [[calculation.level]]
    !> with a topology-independent method (e.g. gfn2, the common choice) and
    !> refine="sp" to re-score the whole final ensemble on one shared
    !> scale -- see search_newnci_metal.f90 and final_tests/
    !> BailarTwist-CoEn3/native_seeded_refined/.
    logical,allocatable :: seed(:)
    !> pi-coordination (eta-bonded) support: an explicit, user-given list
    !> of atom indices forming a coordinating ring/diene face (donor
    !> detection never finds these on its own -- it is hardcoded to the
    !> usual Lewis-basic heteroatoms, N/O/S/P/F/Cl/Br/I, and explicitly
    !> excludes carbon, see search_newnci_metal.f90). When given, a
    !> centroidrange constraint (constraints.f90) is added between this
    !> spec's resolved metal atom(s) and the LIVE centroid of pi_atoms(k)
    !> %idx, giving eta-bonded ligands (COD, Cp, arenes, ...) the same
    !> "survive the MTD bias potential" dissociation protection normal
    !> donors get -- distance-only for now, no angle restraint. Unlike
    !> every other field here, this one is intentionally NOT resolved
    !> through the atom-indexed "last spec wins" mechanism (see
    !> crest_search_newnci_metal): a metal with more than one separate
    !> pi-ligand (e.g. nickelocene's two independent Cp rings) needs TWO
    !> separate [[calculation.nci_metal]] blocks, same atom=, each its
    !> own distinct pi_atoms -- each becomes its own constraint, not
    !> merged into one wrong combined centroid.
    type(pi_atomlist),allocatable :: pi_atoms(:)
  contains
    procedure :: add_spec => nci_metal_settings_add_spec
  end type nci_metal_settings

!========================================================================================!
contains
!========================================================================================!

!========================================================================================!
!> Record one (atom spec, UFF4MOF label) pair, as read from a
!> [[calculation.nci_metal]] block. Specs are not deduplicated: if
!> two specs overlap, resolution order (see crest_search_newnci_metal)
!> decides which one applies to a given atom.
!========================================================================================!
  subroutine nci_metal_settings_add_spec(self,spec,label,restrain,auto,force,angle_fc,r1,theta0,seed,pi_atoms_in)
    class(nci_metal_settings),intent(inout) :: self
    character(len=*),intent(in) :: spec
    character(len=*),intent(in) :: label
    logical,intent(in),optional :: restrain
    logical,intent(in),optional :: auto
    logical,intent(in),optional :: force
    real(wp),intent(in),optional :: angle_fc
    real(wp),intent(in),optional :: r1
    real(wp),intent(in),optional :: theta0
    logical,intent(in),optional :: seed
    integer,intent(in),optional :: pi_atoms_in(:)
    character(len=64),allocatable :: newspec(:)
    character(len=8),allocatable :: newuff(:)
    logical,allocatable :: newrestrain(:),newauto(:),newforce(:),newseed(:)
    real(wp),allocatable :: newanglefc(:),newr1(:),newtheta0(:)
    type(pi_atomlist),allocatable :: newpiatoms(:)
    integer :: ii
    logical :: restrainval,autoval,forceval,seedval
    real(wp) :: anglefcval,r1val,theta0val

    restrainval = .true.
    if (present(restrain)) restrainval = restrain
    autoval = .false.
    if (present(auto)) autoval = auto
    forceval = .false.
    if (present(force)) forceval = force
    anglefcval = -1.0_wp !> <=0 signals "not set, use constraints.f90's own default"
    if (present(angle_fc)) anglefcval = angle_fc
    r1val = -1.0_wp !> <=0 signals "not set, use normal UFF4MOF resolution"
    if (present(r1)) r1val = r1
    theta0val = -1.0_wp !> <=0 signals "not set, use normal UFF4MOF resolution"
    if (present(theta0)) theta0val = theta0
!>--- unlike restrain above, seed defaults to false in every case (even
!>    inside an explicit [[calculation.nci_metal]] block) -- it's the
!>    isomer-generation + enantiomer-duplicate-filter machinery, an
!>    opt-in feature, not a side effect of merely configuring a metal.
    seedval = .false.
    if (present(seed)) seedval = seed

    allocate (newspec(self%nassign+1))
    allocate (newuff(self%nassign+1))
    allocate (newrestrain(self%nassign+1))
    allocate (newauto(self%nassign+1))
    allocate (newforce(self%nassign+1))
    allocate (newanglefc(self%nassign+1))
    allocate (newr1(self%nassign+1))
    allocate (newtheta0(self%nassign+1))
    allocate (newseed(self%nassign+1))
    allocate (newpiatoms(self%nassign+1))
    if (self%nassign > 0) then
      newspec(1:self%nassign) = self%spec(1:self%nassign)
      newuff(1:self%nassign) = self%uff(1:self%nassign)
      newrestrain(1:self%nassign) = self%restrain(1:self%nassign)
      newauto(1:self%nassign) = self%auto(1:self%nassign)
      newforce(1:self%nassign) = self%force(1:self%nassign)
      newanglefc(1:self%nassign) = self%angle_fc(1:self%nassign)
      newr1(1:self%nassign) = self%r1(1:self%nassign)
      newtheta0(1:self%nassign) = self%theta0(1:self%nassign)
      newseed(1:self%nassign) = self%seed(1:self%nassign)
      if (allocated(self%pi_atoms)) then
        do ii = 1,self%nassign
          if (allocated(self%pi_atoms(ii)%idx)) newpiatoms(ii)%idx = self%pi_atoms(ii)%idx
        end do
      end if
    end if
    newspec(self%nassign+1) = adjustl(spec)
    newuff(self%nassign+1) = adjustl(label)
    newrestrain(self%nassign+1) = restrainval
    newauto(self%nassign+1) = autoval
    newforce(self%nassign+1) = forceval
    newanglefc(self%nassign+1) = anglefcval
    newr1(self%nassign+1) = r1val
    newtheta0(self%nassign+1) = theta0val
    newseed(self%nassign+1) = seedval
    if (present(pi_atoms_in)) newpiatoms(self%nassign+1)%idx = pi_atoms_in
    call move_alloc(newspec,self%spec)
    call move_alloc(newuff,self%uff)
    call move_alloc(newanglefc,self%angle_fc)
    call move_alloc(newr1,self%r1)
    call move_alloc(newtheta0,self%theta0)
    call move_alloc(newforce,self%force)
    call move_alloc(newrestrain,self%restrain)
    call move_alloc(newauto,self%auto)
    call move_alloc(newseed,self%seed)
    call move_alloc(newpiatoms,self%pi_atoms)
    self%nassign = self%nassign+1
  end subroutine nci_metal_settings_add_spec

!========================================================================================!
end module nci_metal_settings_mod
!========================================================================================!
!========================================================================================!
