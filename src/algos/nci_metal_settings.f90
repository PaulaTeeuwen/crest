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
  subroutine nci_metal_settings_add_spec(self,spec,label,restrain,auto,force,angle_fc,r1,theta0)
    class(nci_metal_settings),intent(inout) :: self
    character(len=*),intent(in) :: spec
    character(len=*),intent(in) :: label
    logical,intent(in),optional :: restrain
    logical,intent(in),optional :: auto
    logical,intent(in),optional :: force
    real(wp),intent(in),optional :: angle_fc
    real(wp),intent(in),optional :: r1
    real(wp),intent(in),optional :: theta0
    character(len=64),allocatable :: newspec(:)
    character(len=8),allocatable :: newuff(:)
    logical,allocatable :: newrestrain(:),newauto(:),newforce(:)
    real(wp),allocatable :: newanglefc(:),newr1(:),newtheta0(:)
    logical :: restrainval,autoval,forceval
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

    allocate (newspec(self%nassign+1))
    allocate (newuff(self%nassign+1))
    allocate (newrestrain(self%nassign+1))
    allocate (newauto(self%nassign+1))
    allocate (newforce(self%nassign+1))
    allocate (newanglefc(self%nassign+1))
    allocate (newr1(self%nassign+1))
    allocate (newtheta0(self%nassign+1))
    if (self%nassign > 0) then
      newspec(1:self%nassign) = self%spec(1:self%nassign)
      newuff(1:self%nassign) = self%uff(1:self%nassign)
      newrestrain(1:self%nassign) = self%restrain(1:self%nassign)
      newauto(1:self%nassign) = self%auto(1:self%nassign)
      newforce(1:self%nassign) = self%force(1:self%nassign)
      newanglefc(1:self%nassign) = self%angle_fc(1:self%nassign)
      newr1(1:self%nassign) = self%r1(1:self%nassign)
      newtheta0(1:self%nassign) = self%theta0(1:self%nassign)
    end if
    newspec(self%nassign+1) = adjustl(spec)
    newuff(self%nassign+1) = adjustl(label)
    newrestrain(self%nassign+1) = restrainval
    newauto(self%nassign+1) = autoval
    newforce(self%nassign+1) = forceval
    newanglefc(self%nassign+1) = anglefcval
    newr1(self%nassign+1) = r1val
    newtheta0(self%nassign+1) = theta0val
    call move_alloc(newspec,self%spec)
    call move_alloc(newuff,self%uff)
    call move_alloc(newanglefc,self%angle_fc)
    call move_alloc(newr1,self%r1)
    call move_alloc(newtheta0,self%theta0)
    call move_alloc(newforce,self%force)
    call move_alloc(newrestrain,self%restrain)
    call move_alloc(newauto,self%auto)
    self%nassign = self%nassign+1
  end subroutine nci_metal_settings_add_spec

!========================================================================================!
end module nci_metal_settings_mod
!========================================================================================!
!========================================================================================!
