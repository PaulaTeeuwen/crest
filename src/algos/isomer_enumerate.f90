!================================================================================!
! This file is part of crest.
!
! crest is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!================================================================================!
!> Fortran port of tools/isomer_enumerate.py's geometric core: enumerate
!> the symmetry-distinct ways to place a multiset of ligands onto a fixed
!> coordination geometry's sites, so each distinct isomer can seed its
!> own independent nci_metal search instead of relying on one
!> unconstrained run's MTD+GC sampling to discover all of them from a
!> single starting point.
!>
!> Two site-assignments are the SAME isomer iff some PROPER rotation of
!> the coordination polyhedron maps one onto the other. Reflections are
!> deliberately excluded (a proper rotation preserves chirality, so a
!> Bailar-twist Delta/Lambda pair stays two *distinct* isomers). The
!> rotation group is found numerically: brute-force over all site
!> permutations, keeping only those realized by an actual rigid rotation
!> -- using crest's own existing quaternion-based optimal-rotation fit
!> (sorting/ls_rmsd.f90's rmsd()), which can only ever return a PROPER
!> rotation by construction (a quaternion cannot represent a reflection),
!> so no separate determinant check is needed here (the Python version,
!> built on a generic SVD Kabsch fit, needed one explicitly).
!================================================================================!
module isomer_enumerate_mod
  use crest_parameters,only:wp
  use ls_rmsd,only:rmsd
  implicit none
  private

  public :: get_geometry_sites
  public :: find_rotation_group
  public :: n_geometries,geometry_name
  public :: chelate_compatible
  public :: distinct_isomers
  public :: canonical_form

  integer,parameter :: max_raw = 5000   !> generous cap on raw (pre-dedup)
                                         !> assignments; every system this
                                         !> tool has been used on this
                                         !> project stays well under 100

  integer,parameter :: n_geometries = 5
  character(len=24),parameter :: geometry_name(n_geometries) = [ &
    & character(len=24) :: 'tetrahedral','square_planar','trigonal_bipyramidal', &
    & 'square_pyramidal','octahedral']

contains

!========================================================================!
!> Returns the (3,nsites) unit site-vectors for a named geometry, metal
!> at the origin. ok=.false. for an unrecognized name.
  subroutine get_geometry_sites(name,sites,nsites,ok)
    character(len=*),intent(in) :: name
    real(wp),allocatable,intent(out) :: sites(:,:)
    integer,intent(out) :: nsites
    logical,intent(out) :: ok
    real(wp) :: s3
    integer :: i
    ok = .true.
    s3 = sqrt(3.0_wp)/2.0_wp
    select case (trim(name))
    case ('tetrahedral')
      nsites = 4
      allocate (sites(3,nsites))
      sites(:,1) = [1.0_wp,1.0_wp,1.0_wp]
      sites(:,2) = [1.0_wp,-1.0_wp,-1.0_wp]
      sites(:,3) = [-1.0_wp,1.0_wp,-1.0_wp]
      sites(:,4) = [-1.0_wp,-1.0_wp,1.0_wp]
    case ('square_planar')
      nsites = 4
      allocate (sites(3,nsites))
      sites(:,1) = [1.0_wp,0.0_wp,0.0_wp]
      sites(:,2) = [0.0_wp,1.0_wp,0.0_wp]
      sites(:,3) = [-1.0_wp,0.0_wp,0.0_wp]
      sites(:,4) = [0.0_wp,-1.0_wp,0.0_wp]
    case ('trigonal_bipyramidal')
      nsites = 5
      allocate (sites(3,nsites))
      sites(:,1) = [0.0_wp,0.0_wp,1.0_wp]
      sites(:,2) = [0.0_wp,0.0_wp,-1.0_wp]
      sites(:,3) = [1.0_wp,0.0_wp,0.0_wp]
      sites(:,4) = [-0.5_wp,s3,0.0_wp]
      sites(:,5) = [-0.5_wp,-s3,0.0_wp]
    case ('square_pyramidal')
      nsites = 5
      allocate (sites(3,nsites))
      sites(:,1) = [0.0_wp,0.0_wp,1.0_wp]
      sites(:,2) = [1.0_wp,0.0_wp,0.0_wp]
      sites(:,3) = [0.0_wp,1.0_wp,0.0_wp]
      sites(:,4) = [-1.0_wp,0.0_wp,0.0_wp]
      sites(:,5) = [0.0_wp,-1.0_wp,0.0_wp]
    case ('octahedral')
      nsites = 6
      allocate (sites(3,nsites))
      sites(:,1) = [1.0_wp,0.0_wp,0.0_wp]
      sites(:,2) = [-1.0_wp,0.0_wp,0.0_wp]
      sites(:,3) = [0.0_wp,1.0_wp,0.0_wp]
      sites(:,4) = [0.0_wp,-1.0_wp,0.0_wp]
      sites(:,5) = [0.0_wp,0.0_wp,1.0_wp]
      sites(:,6) = [0.0_wp,0.0_wp,-1.0_wp]
    case default
      ok = .false.
      nsites = 0
      allocate (sites(3,0))
      return
    end select
    do i = 1,nsites
      sites(:,i) = sites(:,i)/norm2(sites(:,i))
    end do
  end subroutine get_geometry_sites

!========================================================================!
!> The proper-rotation symmetry group of `sites`, as the set of site-
!> index permutations it induces. Brute-force over all n! permutations
!> (n<=6 here, so <=720 candidates): for each, fit the optimal rotation
!> mapping sites -> sites(perm) via crest's own rmsd() (quaternion-based,
!> always proper), keep the permutation only if that rotation reproduces
!> it essentially exactly (near-zero fit error) -- i.e. the permutation
!> really is realized by a rigid rotation, not just some arbitrary
!> relabeling.
!> `group` is returned as a (nsites,ngroup) integer array, one
!> permutation per column, 1-based site indices.
  subroutine find_rotation_group(sites,nsites,group,ngroup,tol)
    real(wp),intent(in) :: sites(3,nsites)
    integer,intent(in) :: nsites
    integer,allocatable,intent(out) :: group(:,:)
    integer,intent(out) :: ngroup
    real(wp),intent(in),optional :: tol
    real(wp) :: thr
    integer :: perm(nsites),tmp(nsites,720)
    integer :: i,k
    logical :: more
    real(wp) :: Q(3,nsites),U(3,3),xc(3),yc(3),err,g(3,nsites)

    thr = 1.0e-3_wp
    if (present(tol)) thr = tol

    ngroup = 0
    do i = 1,nsites
      perm(i) = i
    end do
    more = .true.
    do while (more)
      do k = 1,nsites
        Q(:,k) = sites(:,perm(k))
      end do
      call rmsd(nsites,sites,Q,1,U,xc,yc,err,.false.,g)
      !> err is the LS-fit RMSD between sites (rotated) and Q under the
      !> trial correspondence given by `perm` -- near-zero iff this
      !> permutation is genuinely realized by a rigid rotation of the
      !> point set, not an arbitrary relabeling that merely happens to
      !> conserve the multiset of site vectors.
      if (err < thr) then
        ngroup = ngroup+1
        tmp(:,ngroup) = perm
      end if
      call next_permutation(perm,nsites,more)
    end do
    allocate (group(nsites,ngroup))
    group = tmp(:,1:ngroup)
  end subroutine find_rotation_group

!========================================================================!
!> Standard lexicographic next-permutation (Narayana's algorithm).
!> more=.false. once `perm` was already the last (descending) one.
  subroutine next_permutation(perm,n,more)
    integer,intent(inout) :: perm(n)
    integer,intent(in) :: n
    logical,intent(out) :: more
    integer :: i,j,t
    i = n-1
    do while (i >= 1)
      if (perm(i) < perm(i+1)) exit
      i = i-1
    end do
    if (i < 1) then
      more = .false.
      return
    end if
    j = n
    do while (perm(j) <= perm(i))
      j = j-1
    end do
    t = perm(i); perm(i) = perm(j); perm(j) = t
    !> reverse the tail after i (length n-i, so (n-i)/2 swaps needed --
    !> NOT (n-i-1)/2, which undercounts by one swap whenever the tail
    !> length is even, silently leaving it half-reversed)
    do j = 0,(n-i)/2-1
      t = perm(i+1+j)
      perm(i+1+j) = perm(n-j)
      perm(n-j) = t
    end do
    more = .true.
  end subroutine next_permutation

!========================================================================!
!> Site-index pairs (1-based) whose vertex angle (metal-site_i /
!> metal-site_j) is <= max_bite_deg -- close enough together to
!> plausibly be bridged by one bidentate ligand.
  subroutine chelate_compatible(sites,nsites,max_bite_deg,compat)
    real(wp),intent(in) :: sites(3,nsites)
    integer,intent(in) :: nsites
    real(wp),intent(in) :: max_bite_deg
    logical,intent(out) :: compat(nsites,nsites)
    integer :: i,j
    real(wp) :: ang,rad2deg
    rad2deg = 180.0_wp/acos(-1.0_wp)
    compat = .false.
    do i = 1,nsites
      compat(i,i) = .true.
      do j = i+1,nsites
        ang = acos(max(-1.0_wp,min(1.0_wp,dot_product(sites(:,i),sites(:,j)))))*rad2deg
        if (ang <= max_bite_deg) then
          compat(i,j) = .true.
          compat(j,i) = .true.
        end if
      end do
    end do
  end subroutine chelate_compatible

!========================================================================!
!> All ways to partition nsites sites into groups matching the ligand
!> units (unit_label(k), unit_denticity(k) for k=1..n_units -- e.g. 3
!> identical bidentate "en" units all share one label code). A
!> denticity-1 unit takes any free site; a denticity-2 unit can only
!> take a free PAIR that is chelate-compatible. For a denticity-2 unit
!> whose label is asymmetric (unit_asym(k)), BOTH site orderings are
!> emitted as separate raw candidates (which end sits where matters);
!> for a symmetric one, only one (order is not physically meaningful).
!> Raw results are returned packed into raw(1:3,1:nsites,1:nraw):
!>   row 1 = the label occupying each site
!>   row 2 = PARTNER site index for a denticity-2 unit's other donor (0
!>           for denticity-1) -- this is NOT redundant with row 1 even
!>           when every unit shares one label (e.g. 3 identical "en"
!>           ligands): row 1 alone is then constant everywhere and
!>           carries no information at all about which two sites are
!>           chelate-bonded TOGETHER as one physical ring, which is
!>           exactly what determines Bailar-twist-type chirality. Row 2
!>           is what restores that pairing/connectivity information
!>           (lost in the Python version's set-of-tuples representation
!>           only because a tuple keeps its two site indices grouped
!>           together already; this flat per-site array needs an
!>           explicit pointer to the same effect).
!>   row 3 = donor END marker (1 or 2) for an asymmetric denticity-2
!>           unit's two sites (which specific donor atom, e.g. glycinate's
!>           N vs O, sits where); irrelevant/always 1 for monodentate or
!>           symmetric denticity-2 units.
  subroutine enumerate_assignments(nsites,n_units,unit_label,unit_denticity,unit_asym, &
                                    compat,raw,nraw)
    integer,intent(in) :: nsites,n_units
    integer,intent(in) :: unit_label(n_units),unit_denticity(n_units)
    logical,intent(in) :: unit_asym(n_units)
    logical,intent(in) :: compat(nsites,nsites)
    integer,intent(out) :: raw(3,nsites,max_raw)
    integer,intent(out) :: nraw
    logical :: used(nsites)
    integer :: cur(3,nsites)
    nraw = 0
    used = .false.
    cur = 0
    call backtrack(1)
  contains
    recursive subroutine backtrack(unit_idx)
      integer,intent(in) :: unit_idx
      integer :: dent,lbl,i,j,a,b,nfree,free_sites(nsites)
      if (unit_idx > n_units) then
        nraw = nraw+1
        if (nraw <= max_raw) raw(:,:,nraw) = cur
        return
      end if
      lbl = unit_label(unit_idx)
      dent = unit_denticity(unit_idx)
      nfree = 0
      do i = 1,nsites
        if (.not.used(i)) then
          nfree = nfree+1
          free_sites(nfree) = i
        end if
      end do
      if (dent == 1) then
        do i = 1,nfree
          a = free_sites(i)
          used(a) = .true.
          cur(1,a) = lbl; cur(2,a) = 0; cur(3,a) = 1
          call backtrack(unit_idx+1)
          used(a) = .false.
        end do
      else
        !> denticity 2 (the only case beyond 1 this project needs)
        do i = 1,nfree
          do j = i+1,nfree
            a = free_sites(i); b = free_sites(j)
            if (.not.compat(a,b)) cycle
            used(a) = .true.; used(b) = .true.
            cur(1,a) = lbl; cur(1,b) = lbl
            cur(2,a) = b; cur(2,b) = a   !> partner pointers, both directions
            !> ordering 1
            cur(3,a) = 1; cur(3,b) = merge(2,1,unit_asym(unit_idx))
            call backtrack(unit_idx+1)
            if (unit_asym(unit_idx)) then
              !> ordering 2: genuinely different for an asymmetric
              !> ligand (which end sits at a vs b); a no-op re-visit for
              !> a symmetric one, correctly collapsed later
              cur(3,a) = 2; cur(3,b) = 1
              call backtrack(unit_idx+1)
            end if
            used(a) = .false.; used(b) = .false.
          end do
        end do
      end if
    end subroutine backtrack
  end subroutine enumerate_assignments

!========================================================================!
!> The lexicographically-smallest image of one (label,partner,end)
!> assignment under the rotation group -- two raw assignments are the
!> SAME distinct isomer iff they share a canonical form. Row 1 (label)
!> is compared first, then row 2 (partner site, REMAPPED through the
!> same rotation -- not copied raw, since a site index that was
!> meaningful before the rotation must stay meaningful after it: if site
!> s's partner was site p, then after rotation g the ligand now at
!> site group(s,g) has its partner at site group(p,g), not at the old
!> numeric value p), then row 3 (donor end) as the final tie-break.
  subroutine canonical_form(assign,nsites,group,ngroup,canon)
    integer,intent(in) :: nsites,ngroup
    integer,intent(in) :: assign(3,nsites)
    integer,intent(in) :: group(nsites,ngroup)
    integer,intent(out) :: canon(3,nsites)
    integer :: g,s,trial(3,nsites)
    logical :: better,tie
    canon = 0
    do g = 1,ngroup
      do s = 1,nsites
        trial(1,group(s,g)) = assign(1,s)
        trial(3,group(s,g)) = assign(3,s)
        if (assign(2,s) > 0) then
          trial(2,group(s,g)) = group(assign(2,s),g)
        else
          trial(2,group(s,g)) = 0
        end if
      end do
      if (g == 1) then
        canon = trial
      else
        better = .false.
        tie = .true.
        do s = 1,nsites
          if (trial(1,s) < canon(1,s)) then
            better = .true.; tie = .false.; exit
          else if (trial(1,s) > canon(1,s)) then
            tie = .false.; exit
          end if
        end do
        if (tie) then
          do s = 1,nsites
            if (trial(2,s) < canon(2,s)) then
              better = .true.; tie = .false.; exit
            else if (trial(2,s) > canon(2,s)) then
              tie = .false.; exit
            end if
          end do
        end if
        if (tie) then
          do s = 1,nsites
            if (trial(3,s) < canon(3,s)) then
              better = .true.; exit
            else if (trial(3,s) > canon(3,s)) then
              exit
            end if
          end do
        end if
        if (better) canon = trial
      end if
    end do
  end subroutine canonical_form

!========================================================================!
!> Top-level entry point: every symmetry-distinct isomer for placing
!> n_units ligands (unit_label/unit_denticity/unit_asym, as in
!> enumerate_assignments) onto `geomname`'s coordination sites.
!> Returns distinct canonical assignments packed into
!> isomers(1:3,1:nsites,1:niso), deduplicated and sorted.
  subroutine distinct_isomers(geomname,n_units,unit_label,unit_denticity,unit_asym, &
                               max_bite_deg,isomers,niso,ok)
    character(len=*),intent(in) :: geomname
    integer,intent(in) :: n_units
    integer,intent(in) :: unit_label(n_units),unit_denticity(n_units)
    logical,intent(in) :: unit_asym(n_units)
    real(wp),intent(in) :: max_bite_deg
    integer,allocatable,intent(out) :: isomers(:,:,:)
    integer,intent(out) :: niso
    logical,intent(out) :: ok
    real(wp),allocatable :: sites(:,:)
    integer :: nsites
    integer,allocatable :: group(:,:)
    integer :: ngroup
    logical,allocatable :: compat(:,:)
    integer,allocatable :: raw(:,:,:)
    integer :: nraw
    integer,allocatable :: canon_all(:,:,:)
    integer :: i,k,idx(max_raw)
    logical :: isdup

    call get_geometry_sites(geomname,sites,nsites,ok)
    if (.not.ok) return
    if (sum(unit_denticity) /= nsites) then
      ok = .false.
      niso = 0
      allocate (isomers(3,nsites,0))
      return
    end if
    call find_rotation_group(sites,nsites,group,ngroup)
    allocate (compat(nsites,nsites))
    call chelate_compatible(sites,nsites,max_bite_deg,compat)

    allocate (raw(3,nsites,max_raw))
    call enumerate_assignments(nsites,n_units,unit_label,unit_denticity,unit_asym, &
                                compat,raw,nraw)
    if (nraw > max_raw) nraw = max_raw  !> safety; should never trigger in practice

    allocate (canon_all(3,nsites,nraw))
    do i = 1,nraw
      call canonical_form(raw(:,:,i),nsites,group,ngroup,canon_all(:,:,i))
    end do

    !> dedup: keep the first occurrence of each distinct canonical form
    !> (O(n^2) over nraw, fine at this scale -- nraw stays in the tens
    !> to low hundreds for every system this tool is used on)
    niso = 0
    do i = 1,nraw
      isdup = .false.
      do k = 1,niso
        if (all(canon_all(:,:,i) == canon_all(:,:,idx(k)))) then
          isdup = .true.
          exit
        end if
      end do
      if (.not.isdup) then
        niso = niso+1
        idx(niso) = i
      end if
    end do

    allocate (isomers(3,nsites,niso))
    do i = 1,niso
      isomers(:,:,i) = canon_all(:,:,idx(i))
    end do
    ok = .true.
  end subroutine distinct_isomers

end module isomer_enumerate_mod
