!================================================================================!
! This file is part of crest.
!
! crest is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!================================================================================!
!> Fortran port of tools/generate_isomer_structures.py: given a real
!> starting structure for a metal complex, detect its ligand fragments
!> and their denticity from connectivity, enumerate the symmetry-distinct
!> ways those SAME fragments could be arranged on the coordination
!> polyhedron (isomer_enumerate_mod), and for every isomer other than the
!> input's own, rigidly transplant each ligand fragment (whole fragment,
!> not just its donor atom) onto its new target site direction(s) --
!> producing one full starting structure per distinct isomer, for a
!> per-isomer seeded nci_metal search.
!>
!> All internal geometry here is done in Angstrom (mirroring the
!> reference Python tool exactly, including its hand-tuned clash/
!> expansion distances), converted from/to crest's Bohr convention only
!> at the module boundary.
!================================================================================!
module isomer_generate_mod
  use crest_parameters,only:wp,stdout,autoaa
  use isomer_enumerate_mod,only:get_geometry_sites,find_rotation_group,canonical_form, &
    & chelate_compatible,distinct_isomers
  use ls_rmsd,only:rmsd
  implicit none
  private

  public :: generate_isomer_structures

  integer,parameter :: max_frag = 30
  integer,parameter :: max_fatoms = 500
  integer,parameter :: max_isomers = 200
  integer,parameter :: max_wl_deg = 12

  type :: fragment_t
    integer :: natoms = 0
    integer :: atoms(max_fatoms) = 0
    integer :: ndonor = 0
    integer :: donors(2) = 0     !> this project never needs denticity > 2
    integer :: denticity = 0
    integer :: label = 0
    logical :: asym = .false.
  end type fragment_t

contains

!========================================================================!
!> Local covalent radii (Angstrom), atomic-number indexed, matching the
!> reference Python tool's COVALENT_RADIUS/METAL_RADIUS tables exactly
!> (kept local/independent of crest's own rcov table in miscdata.f90,
!> which is pre-scaled by 4/3 for its own damping-function purpose and
!> would silently double-scale a hard distance cutoff used here).
  pure function cov_radius(iz) result(r)
    integer,intent(in) :: iz
    real(wp) :: r
    select case (iz)
    case (1); r = 0.31_wp      !> H
    case (6); r = 0.76_wp      !> C
    case (7); r = 0.71_wp      !> N
    case (8); r = 0.66_wp      !> O
    case (9); r = 0.57_wp      !> F
    case (15); r = 1.07_wp     !> P
    case (16); r = 1.05_wp     !> S
    case (17); r = 1.02_wp     !> Cl
    case (35); r = 1.20_wp     !> Br
    case (53); r = 1.39_wp     !> I
    case (24); r = 1.39_wp     !> Cr
    case (25); r = 1.39_wp     !> Mn
    case (26); r = 1.32_wp     !> Fe
    case (27); r = 1.26_wp     !> Co
    case (28); r = 1.24_wp     !> Ni
    case (29); r = 1.32_wp     !> Cu
    case (30); r = 1.22_wp     !> Zn
    case (46); r = 1.39_wp     !> Pd
    case (78); r = 1.36_wp     !> Pt
    case default; r = 0.75_wp
    end select
  end function cov_radius

!========================================================================!
  subroutine detect_bonds(nat,iz,xyz_aa,bonds,scale)
    integer,intent(in) :: nat,iz(nat)
    real(wp),intent(in) :: xyz_aa(3,nat)
    logical,intent(out) :: bonds(nat,nat)
    real(wp),intent(in),optional :: scale
    real(wp) :: sc,cutoff
    integer :: i,j
    sc = 1.2_wp
    if (present(scale)) sc = scale
    bonds = .false.
    do i = 1,nat
      do j = i+1,nat
        cutoff = sc*(cov_radius(iz(i))+cov_radius(iz(j)))
        if (norm2(xyz_aa(:,i)-xyz_aa(:,j)) <= cutoff) then
          bonds(i,j) = .true.; bonds(j,i) = .true.
        end if
      end do
    end do
  end subroutine detect_bonds

!========================================================================!
!> One fragment per ligand: BFS/DFS from each metal-neighbor, excluding
!> the metal itself, over the bond graph.
  subroutine detect_fragments(nat,metal_idx,bonds,fragments,nfrag)
    integer,intent(in) :: nat,metal_idx
    logical,intent(in) :: bonds(nat,nat)
    type(fragment_t),intent(out) :: fragments(max_frag)
    integer,intent(out) :: nfrag
    logical :: assigned(nat),seen(nat)
    integer :: stackarr(nat),nstack
    integer :: start,cur,nxt,i,k,nd
    assigned = .false.
    assigned(metal_idx) = .true.
    nfrag = 0
    do start = 1,nat
      if (start == metal_idx) cycle
      if (.not.bonds(metal_idx,start)) cycle
      if (assigned(start)) cycle
      seen = .false.
      seen(start) = .true.
      nstack = 1
      stackarr(1) = start
      do while (nstack > 0)
        cur = stackarr(nstack)
        nstack = nstack-1
        do nxt = 1,nat
          if (nxt == metal_idx.or.seen(nxt)) cycle
          if (.not.bonds(cur,nxt)) cycle
          seen(nxt) = .true.
          nstack = nstack+1
          stackarr(nstack) = nxt
        end do
      end do
      nfrag = nfrag+1
      fragments(nfrag)%natoms = 0
      do i = 1,nat
        if (seen(i)) then
          fragments(nfrag)%natoms = fragments(nfrag)%natoms+1
          fragments(nfrag)%atoms(fragments(nfrag)%natoms) = i
        end if
      end do
      nd = 0
      do k = 1,fragments(nfrag)%natoms
        i = fragments(nfrag)%atoms(k)
        if (bonds(metal_idx,i)) then
          nd = nd+1
          if (nd <= 2) fragments(nfrag)%donors(nd) = i
        end if
      end do
      fragments(nfrag)%ndonor = nd
      fragments(nfrag)%denticity = nd
      assigned = assigned.or.seen
    end do
  end subroutine detect_fragments

!========================================================================!
!> ascending insertion sort of an integer array's first n entries
  subroutine isort(a,n)
    integer,intent(inout) :: a(:)
    integer,intent(in) :: n
    integer :: i,j,t
    do i = 2,n
      t = a(i); j = i-1
      do while (j >= 1)
        if (a(j) <= t) exit
        a(j+1) = a(j); j = j-1
      end do
      a(j+1) = t
    end do
  end subroutine isort

!========================================================================!
!> 1-WL-style color refinement over the bonding graph restricted to
!> `atoms` (element labels + connectivity only). Two atoms end up with
!> the same final color only if no local structural feature (out to
!> natoms rounds away) distinguishes them -- used both to type ligand
!> fragments by bonding topology (not just elemental formula, which
!> would wrongly conflate two constitutionally different ligands that
!> happen to share an atom count) and to detect whether a bidentate
!> fragment's two donors are chemically interchangeable.
  subroutine color_refine(natoms,atoms,nat,iz,bonds,color)
    integer,intent(in) :: natoms,atoms(natoms),nat,iz(nat)
    logical,intent(in) :: bonds(nat,nat)
    integer,intent(out) :: color(natoms)
    integer :: key(max_wl_deg+1,natoms),uniq(max_wl_deg+1,max_fatoms)
    integer :: nbr(max_wl_deg)
    integer :: round,i,j,a,b,deg,nuniq,k
    logical :: found
    do i = 1,natoms
      color(i) = iz(atoms(i))
    end do
    do round = 1,natoms
      do i = 1,natoms
        a = atoms(i)
        deg = 0
        do j = 1,natoms
          if (i == j) cycle
          b = atoms(j)
          if (bonds(a,b)) then
            deg = deg+1
            if (deg <= max_wl_deg) nbr(deg) = color(j)
          end if
        end do
        if (deg > max_wl_deg) deg = max_wl_deg  !> never expected for real ligands
        call isort(nbr,deg)
        key(1,i) = color(i)
        do k = 1,max_wl_deg
          if (k <= deg) then
            key(1+k,i) = nbr(k)
          else
            key(1+k,i) = -1
          end if
        end do
      end do
      nuniq = 0
      do i = 1,natoms
        found = .false.
        do k = 1,nuniq
          if (all(key(:,i) == uniq(:,k))) then
            found = .true.; exit
          end if
        end do
        if (.not.found) then
          nuniq = nuniq+1
          uniq(:,nuniq) = key(:,i)
        end if
      end do
      !> sort the unique keys lexicographically so the rank assigned
      !> below is a deterministic canonical color, not iteration-order
      !> dependent
      call sort_key_columns(uniq,nuniq,max_wl_deg+1)
      do i = 1,natoms
        do k = 1,nuniq
          if (all(key(:,i) == uniq(:,k))) then
            color(i) = k
            exit
          end if
        end do
      end do
    end do
  end subroutine color_refine

!========================================================================!
!> ascending lexicographic insertion sort of the first ncol columns of
!> an (nrow,*) integer array
  subroutine sort_key_columns(a,ncol,nrow)
    integer,intent(inout) :: a(:,:)
    integer,intent(in) :: ncol,nrow
    integer :: i,j,t(nrow)
    logical :: less
    do i = 2,ncol
      t = a(:,i); j = i-1
      do while (j >= 1)
        less = .false.
        block
          integer :: r
          do r = 1,nrow
            if (t(r) < a(r,j)) then
              less = .true.; exit
            else if (t(r) > a(r,j)) then
              exit
            end if
          end do
        end block
        if (.not.less) exit
        a(:,j+1) = a(:,j); j = j-1
      end do
      a(:,j+1) = t
    end do
  end subroutine sort_key_columns

!========================================================================!
!> Topology-type key for a fragment: sorted (atomic_number,count) list
!> (the formula) followed by the sorted final color-refinement colors.
!> Two fragments are the SAME ligand type iff these keys match exactly.
!> Packed into a single padded integer array for easy comparison.
  subroutine fragment_topokey(frag,nat,iz,bonds,key,keylen)
    type(fragment_t),intent(in) :: frag
    integer,intent(in) :: nat,iz(nat)
    logical,intent(in) :: bonds(nat,nat)
    integer,intent(out) :: key(2*118+max_fatoms)
    integer,intent(out) :: keylen
    integer :: counts(118),z,i,n
    integer :: color(frag%natoms),colorsorted(frag%natoms)
    counts = 0
    do i = 1,frag%natoms
      counts(iz(frag%atoms(i))) = counts(iz(frag%atoms(i)))+1
    end do
    n = 0
    key = 0
    do z = 1,118
      if (counts(z) > 0) then
        key(n+1) = z; key(n+2) = counts(z)
        n = n+2
      end if
    end do
    call color_refine(frag%natoms,frag%atoms,nat,iz,bonds,color)
    colorsorted = color
    call isort(colorsorted,frag%natoms)
    key(n+1:n+frag%natoms) = colorsorted(1:frag%natoms)
    n = n+frag%natoms
    keylen = n
  end subroutine fragment_topokey

!========================================================================!
!> Rodrigues' rotation formula: rotation by `angle` (rad) about `axis`.
  function rotmat_axis_angle(axis,angle) result(R)
    real(wp),intent(in) :: axis(3),angle
    real(wp) :: R(3,3),ax(3),K(3,3),I3(3,3)
    ax = axis/norm2(axis)
    K = reshape([0.0_wp,ax(3),-ax(2), -ax(3),0.0_wp,ax(1), ax(2),-ax(1),0.0_wp],[3,3])
    !> column-major reshape above already yields the standard
    !> skew-symmetric cross-product matrix for ax in Fortran's
    !> column-major storage (K(i,j) as built matches [[0,-az,ay],
    !> [az,0,-ax],[-ay,ax,0]] when read row-wise)
    I3 = 0.0_wp
    I3(1,1) = 1.0_wp; I3(2,2) = 1.0_wp; I3(3,3) = 1.0_wp
    R = I3+sin(angle)*K+(1.0_wp-cos(angle))*matmul(K,K)
  end function rotmat_axis_angle

!========================================================================!
!> The minimal proper rotation mapping unit vector a exactly onto unit
!> vector b (handles the antiparallel case explicitly, matching the
!> Python tool's align_vectors exactly: a 180-deg rotation must be
!> about an axis PERPENDICULAR to a, not about a itself).
  function align_vectors(a,b) result(R)
    real(wp),intent(in) :: a(3),b(3)
    real(wp) :: R(3,3),axis(3),angle,perp(3),I3(3,3)
    axis = cross3(a,b)
    if (norm2(axis) < 1.0e-8_wp) then
      if (dot_product(a,b) > 0.0_wp) then
        R = 0.0_wp; R(1,1) = 1.0_wp; R(2,2) = 1.0_wp; R(3,3) = 1.0_wp
        return
      end if
      if (abs(a(1)) < 0.9_wp) then
        perp = [1.0_wp,0.0_wp,0.0_wp]
      else
        perp = [0.0_wp,1.0_wp,0.0_wp]
      end if
      perp = perp-a*dot_product(perp,a)
      perp = perp/norm2(perp)
      I3 = 0.0_wp; I3(1,1) = 1.0_wp; I3(2,2) = 1.0_wp; I3(3,3) = 1.0_wp
      R = 2.0_wp*outer3(perp,perp)-I3
      return
    end if
    angle = acos(max(-1.0_wp,min(1.0_wp,dot_product(a,b))))
    R = rotmat_axis_angle(axis,angle)
  end function align_vectors

  pure function cross3(a,b) result(c)
    real(wp),intent(in) :: a(3),b(3)
    real(wp) :: c(3)
    c(1) = a(2)*b(3)-a(3)*b(2)
    c(2) = a(3)*b(1)-a(1)*b(3)
    c(3) = a(1)*b(2)-a(2)*b(1)
  end function cross3

  pure function outer3(a,b) result(M)
    real(wp),intent(in) :: a(3),b(3)
    real(wp) :: M(3,3)
    integer :: i,j
    do j = 1,3
      do i = 1,3
        M(i,j) = a(i)*b(j)
      end do
    end do
  end function outer3

!========================================================================!
!> The proper rotation mapping unit vectors a1->b1 AND a2->b2 exactly
!> (ported directly from the Python tool's exact_two_point_rotation --
!> see its docstring there for why this closed-form construction is
!> used instead of an SVD/Kabsch fit).
  function exact_two_point_rotation(a1,b1,a2,b2) result(R)
    real(wp),intent(in) :: a1(3),b1(3),a2(3),b2(3)
    real(wp) :: R(3,3),R1(3,3),R2(3,3),a2p(3),p1(3),p2(3),n1,n2
    R1 = align_vectors(a1,b1)
    a2p = matmul(R1,a2)
    p1 = a2p-b1*dot_product(a2p,b1)
    p2 = b2-b1*dot_product(b2,b1)
    n1 = norm2(p1); n2 = norm2(p2)
    if (n1 < 1.0e-8_wp.or.n2 < 1.0e-8_wp) then
      R = R1
      return
    end if
    p1 = p1/n1; p2 = p2/n2
    if (norm2(cross3(p1,p2)) < 1.0e-8_wp.and.dot_product(p1,p2) < 0.0_wp) then
      R2 = rotmat_axis_angle(b1,acos(-1.0_wp))
    else
      R2 = align_vectors(p1,p2)
    end if
    R = matmul(R2,R1)
  end function exact_two_point_rotation

!========================================================================!
!> Number of (group,other) atom pairs closer than `slack` times their
!> covalent-radius sum.
  function clash_count(nat,iz,xyz_aa,ngroup,group_atoms,nother,other_atoms,slack) result(n)
    integer,intent(in) :: nat,iz(nat),ngroup,group_atoms(ngroup),nother,other_atoms(nother)
    real(wp),intent(in) :: xyz_aa(3,nat)
    real(wp),intent(in),optional :: slack
    integer :: n,i,j,a,b
    real(wp) :: sl,cutoff
    sl = 0.75_wp
    if (present(slack)) sl = slack
    n = 0
    do i = 1,ngroup
      a = group_atoms(i)
      do j = 1,nother
        b = other_atoms(j)
        cutoff = sl*(cov_radius(iz(a))+cov_radius(iz(b)))
        if (norm2(xyz_aa(:,a)-xyz_aa(:,b)) < cutoff) n = n+1
      end do
    end do
  end function clash_count

!========================================================================!
!> Which tabulated geometry (same site count) best matches the current
!> donor-direction set, by comparing sorted pairwise angles.
  subroutine best_matching_geometry(ndonor,donor_dirs,geomname)
    integer,intent(in) :: ndonor
    real(wp),intent(in) :: donor_dirs(3,ndonor)
    character(len=*),intent(out) :: geomname
    character(len=24) :: candidates(5)
    integer :: ncand,ic,nsites
    real(wp),allocatable :: sites(:,:)
    real(wp) :: target_ang(ndonor*(ndonor-1)/2),trial_ang(ndonor*(ndonor-1)/2)
    real(wp) :: err,best_err
    logical :: ok
    candidates = [character(len=24) :: 'tetrahedral','square_planar', &
      & 'trigonal_bipyramidal','square_pyramidal','octahedral']
    call pairwise_angles(ndonor,donor_dirs,target_ang)
    best_err = huge(1.0_wp)
    geomname = ''
    do ic = 1,5
      call get_geometry_sites(trim(candidates(ic)),sites,nsites,ok)
      if (ok.and.nsites == ndonor) then
        call pairwise_angles(nsites,sites,trial_ang)
        err = sum((trial_ang-target_ang)**2)
        if (err < best_err) then
          best_err = err
          geomname = trim(candidates(ic))
        end if
      end if
      if (allocated(sites)) deallocate (sites)
    end do
  end subroutine best_matching_geometry

  subroutine pairwise_angles(n,v,out)
    integer,intent(in) :: n
    real(wp),intent(in) :: v(3,n)
    real(wp),intent(out) :: out(n*(n-1)/2)
    integer :: i,j,k
    real(wp) :: c
    k = 0
    do i = 1,n
      do j = i+1,n
        k = k+1
        c = max(-1.0_wp,min(1.0_wp,dot_product(v(:,i),v(:,j))))
        out(k) = acos(c)*180.0_wp/acos(-1.0_wp)
      end do
    end do
    call rsort(out,k)
  end subroutine pairwise_angles

  subroutine rsort(a,n)
    real(wp),intent(inout) :: a(:)
    integer,intent(in) :: n
    integer :: i,j
    real(wp) :: t
    do i = 2,n
      t = a(i); j = i-1
      do while (j >= 1)
        if (a(j) <= t) exit
        a(j+1) = a(j); j = j-1
      end do
      a(j+1) = t
    end do
  end subroutine rsort

!========================================================================!
!> Top-level entry point. `mol_at`/`mol_xyz_bohr` are the plain arrays
!> from a crest coord object (xyz in Bohr, at = atomic numbers).
!> `metal_idx` is 1-based. Returns every newly-generated isomer's
!> coordinates (Bohr) packed into structures(3,nat,nstruct) -- the
!> input's own arrangement is detected and skipped, exactly as the
!> Python tool does.
  subroutine generate_isomer_structures(nat,at,xyz_bohr,metal_idx,max_bite_deg, &
                                         structures,nstruct,geomname,ok)
    integer,intent(in) :: nat,at(nat),metal_idx
    real(wp),intent(in) :: xyz_bohr(3,nat)
    real(wp),intent(in) :: max_bite_deg
    real(wp),allocatable,intent(out) :: structures(:,:,:)
    integer,intent(out) :: nstruct
    character(len=*),intent(out) :: geomname
    logical,intent(out) :: ok

    real(wp) :: xyz_aa(3,nat),metal_pos(3)
    logical :: bonds(nat,nat)
    type(fragment_t) :: fragments(max_frag)
    integer :: nfrag
    integer :: label_key(max_frag,2*118+max_fatoms),label_keylen(max_frag)
    integer :: label_code(max_frag)
    integer :: nlabels
    integer :: donor_frag(2*max_frag),donor_atom(2*max_frag),ndonor_total
    real(wp) :: donor_dirs(3,2*max_frag)
    integer :: nsites
    real(wp),allocatable :: geom_sites(:,:)
    integer :: donor_to_site(2*max_frag)
    logical :: used_site(2*max_frag)
    integer :: unit_label(max_frag),unit_denticity(max_frag)
    logical :: unit_asym(max_frag)
    integer,allocatable :: isomers_out(:,:,:)
    integer :: niso
    integer,allocatable :: group(:,:)
    integer :: ngroup
    integer :: input_raw(3,2*max_frag),input_canon(3,2*max_frag)
    real(wp) :: new_coords(3,nat)
    integer :: i,j,k,f,s,p,site_consumed_up_to(max_frag)
    integer :: frag_of_label_consumed(max_frag)
    logical :: input_matches
    real(wp) :: c

    ok = .true.
    nstruct = 0
    xyz_aa = xyz_bohr*autoaa
    metal_pos = xyz_aa(:,metal_idx)

    call detect_bonds(nat,at,xyz_aa,bonds)
    call detect_fragments(nat,metal_idx,bonds,fragments,nfrag)
    if (nfrag < 1.or.nfrag > max_frag) then
      ok = .false.
      allocate (structures(3,nat,0))
      return
    end if

    !> type each fragment by bonding topology, and detect donor symmetry
    do f = 1,nfrag
      call fragment_topokey(fragments(f),nat,at,bonds,label_key(f,:),label_keylen(f))
      fragments(f)%asym = .false.
      if (fragments(f)%denticity >= 2) then
        block
          integer :: col(fragments(f)%natoms)
          call color_refine(fragments(f)%natoms,fragments(f)%atoms,nat,at,bonds,col)
          block
            integer :: c1,c2,ii
            c1 = 0; c2 = 0
            do ii = 1,fragments(f)%natoms
              if (fragments(f)%atoms(ii) == fragments(f)%donors(1)) c1 = col(ii)
              if (fragments(f)%atoms(ii) == fragments(f)%donors(2)) c2 = col(ii)
            end do
            fragments(f)%asym = (c1 /= c2)
          end block
        end block
      end if
    end do

    !> assign integer label codes: same key -> same code, first-occurrence order
    nlabels = 0
    do f = 1,nfrag
      fragments(f)%label = 0
      do k = 1,nlabels
        if (label_keylen(f) == label_keylen(label_code(k)).and. &
          & all(label_key(f,1:label_keylen(f)) == label_key(label_code(k),1:label_keylen(f)))) then
          fragments(f)%label = k
          exit
        end if
      end do
      if (fragments(f)%label == 0) then
        nlabels = nlabels+1
        label_code(nlabels) = f
        fragments(f)%label = nlabels
      end if
    end do

    !> collect all donor directions, in fragment-detection order
    ndonor_total = 0
    do f = 1,nfrag
      do j = 1,fragments(f)%denticity
        ndonor_total = ndonor_total+1
        donor_frag(ndonor_total) = f
        donor_atom(ndonor_total) = fragments(f)%donors(j)
        donor_dirs(:,ndonor_total) = (xyz_aa(:,fragments(f)%donors(j))-metal_pos)
        donor_dirs(:,ndonor_total) = donor_dirs(:,ndonor_total)/norm2(donor_dirs(:,ndonor_total))
      end do
    end do

    call best_matching_geometry(ndonor_total,donor_dirs(:,1:ndonor_total),geomname)
    if (len_trim(geomname) == 0) then
      ok = .false.
      allocate (structures(3,nat,0))
      return
    end if
    call get_geometry_sites(trim(geomname),geom_sites,nsites,ok)
    if (.not.ok.or.nsites /= ndonor_total) then
      ok = .false.
      allocate (structures(3,nat,0))
      return
    end if

!>--- ROOT-CAUSE FIX: get_geometry_sites returns fixed, generic axes
!>    (+-x/+-y/+-z-type unit vectors for the named geometry's own idealized
!>    template), in a reference frame that has nothing to do with how this
!>    molecule's REAL donors happen to be oriented in space. Every use of
!>    geom_sites below -- the greedy donor-to-site assignment just below,
!>    and (in generate_isomer_structures) the actual target directions
!>    ligands get rotated onto -- compared/placed vectors from two
!>    UNRELATED frames without this fix, which silently scrambles which
!>    sites are truly antipodal to each other in real 3D space. Confirmed
!>    directly on Co(NO2)3(NH3)3: the facial input structure's own real
!>    NO2 donors were mapped onto a template site SET that included an
!>    antipodal pair (i.e. logically "meridional" in site-space) despite
!>    being genuinely facial in real space -- so a "different, meridional"
!>    site-label pattern generated for a seed candidate never actually
!>    corresponded to a real meridional arrangement, and the candidate's
!>    real geometry was fac from the very first frame, not something that
!>    later relaxed away from mer.
!>
!>    Fix: rotate the WHOLE idealized template (every site in geom_sites)
!>    by the one rigid rotation that exactly maps two of the template's
!>    own sites onto two of this molecule's real donor directions,
!>    choosing whichever site PAIR has the pairwise angle closest to the
!>    real angle between those two donors (exact_two_point_rotation,
!>    already used below for denticity-2 placement, handles the parallel/
!>    antiparallel edge cases). After this, geom_sites lives in the
!>    molecule's own frame, so dot products against real donor_dirs (and
!>    later, real atom placement against geom_sites) are geometrically
!>    meaningful.
!>
!>    The anchor PAIR of real donors must itself be well-conditioned:
!>    if the two chosen donors are (anti)parallel, so is whichever
!>    template site pair best matches them, and exact_two_point_rotation
!>    then has NOTHING left to fix the twist around that shared axis (its
!>    second point maps there automatically once the first is aligned) --
!>    every OTHER site silently lands at an arbitrary rotation around
!>    that axis, corrupting the whole donor-to-site assignment even
!>    though the (anti)parallel pair itself still "looks" aligned.
!>    Confirmed directly: this degeneracy fires whenever the first two
!>    donors in detection order happen to be mutually trans (e.g. the two
!>    Cl- of a square-planar Cu(II) center), and silently discards a real
!>    cis/trans isomer that generate_isomer_structures correctly found
!>    before this whole alignment fix existed. Avoid it by choosing,
!>    among ALL pairs of real donors, the one whose mutual angle is
!>    farthest from both 0 and 180 deg (closest to 90 deg) as the anchor.
!>    A single exact 2-point fit is not enough on its own: it reproduces
!>    those two directions exactly and dumps 100% of whatever distortion
!>    the real (never perfectly idealized) molecule has onto the OTHER
!>    sites, which can misassign them -- confirmed directly: it recovered
!>    the wrong (far higher-energy, visibly distorted) "mer" candidate
!>    for Co(NO2)3(NH3)3 instead of the real one. Used here only to get
!>    an initial, non-degenerate donor<->site correspondence; the ACTUAL
!>    rotation applied is a full least-squares (Kabsch) fit over ALL
!>    donors at once (crest's own ls_rmsd, same routine find_rotation_group
!>    already uses), which spreads real-molecule distortion evenly
!>    instead of concentrating it.
    if (ndonor_total >= 2) then
      block
        integer :: ap,aq,bestap,bestaq,di,dj,bestdi,bestdj
        real(wp) :: ang_donor,ang_site,besterr,thiserr,Ralign(3,3),bestcond
        real(wp) :: geom_sites_orig(3,nsites),dots(nsites)
        integer :: rough_to_site(ndonor_total)
        logical :: rough_used(nsites)
        real(wp) :: U(3,3),xc(3),yc(3),fiterr,g(3,ndonor_total)
        real(wp) :: coord1(3,ndonor_total),coord2(3,ndonor_total)
        geom_sites_orig = geom_sites
!>--- step 1: pick a well-conditioned (least collinear) real donor pair
!>    as the anchor, so the initial exact-fit rotation has no degenerate
!>    (underdetermined-twist) direction -- see the note above this block
!>    in an earlier version of this fix.
        bestdi = 1; bestdj = min(2,ndonor_total)
        bestcond = huge(1.0_wp)
        do di = 1,ndonor_total-1
          do dj = di+1,ndonor_total
            thiserr = abs(dot_product(donor_dirs(:,di),donor_dirs(:,dj)))
            if (thiserr < bestcond) then
              bestcond = thiserr
              bestdi = di; bestdj = dj
            end if
          end do
        end do
        ang_donor = acos(max(-1.0_wp,min(1.0_wp, &
          & dot_product(donor_dirs(:,bestdi),donor_dirs(:,bestdj)))))
        besterr = huge(1.0_wp)
        bestap = 1; bestaq = min(2,nsites)
        do ap = 1,nsites
          do aq = 1,nsites
            if (ap == aq) cycle
            ang_site = acos(max(-1.0_wp,min(1.0_wp, &
              & dot_product(geom_sites(:,ap),geom_sites(:,aq)))))
            thiserr = abs(ang_site-ang_donor)
            if (thiserr < besterr) then
              besterr = thiserr
              bestap = ap; bestaq = aq
            end if
          end do
        end do
        Ralign = exact_two_point_rotation(geom_sites(:,bestap),donor_dirs(:,bestdi), &
          & geom_sites(:,bestaq),donor_dirs(:,bestdj))
        do s = 1,nsites
          geom_sites(:,s) = matmul(Ralign,geom_sites(:,s))
        end do
!>--- step 2: with that rough alignment, get a full donor<->site
!>    correspondence (same greedy rule used for the real assignment
!>    below), then refit ALL of it at once with a proper Kabsch fit.
        rough_used = .false.
        do i = 1,ndonor_total
          do s = 1,nsites
            dots(s) = -dot_product(donor_dirs(:,i),geom_sites(:,s))
          end do
          ap = 0
          do s = 1,nsites
            if (.not.rough_used(s)) then
              if (ap == 0.or.dots(s) < dots(ap)) ap = s
            end if
          end do
          rough_to_site(i) = ap
          rough_used(ap) = .true.
        end do
        do i = 1,ndonor_total
          coord1(:,i) = geom_sites_orig(:,rough_to_site(i))
          coord2(:,i) = donor_dirs(:,i)
        end do
        call rmsd(ndonor_total,coord1,coord2,1,U,xc,yc,fiterr,.false.,g)
        do s = 1,nsites
          geom_sites(:,s) = matmul(U,geom_sites_orig(:,s))
        end do
      end block
    end if

    !> greedy nearest-available-site assignment, same order/criterion as
    !> the Python tool (descending dot product, first free site wins)
    used_site = .false.
    do i = 1,ndonor_total
      block
        integer :: order(nsites),best_s
        real(wp) :: dots(nsites)
        do s = 1,nsites
          dots(s) = -dot_product(donor_dirs(:,i),geom_sites(:,s))
          order(s) = s
        end do
        !> simple selection: pick the unused site with smallest dots (= largest dot product)
        best_s = 0
        do s = 1,nsites
          if (.not.used_site(s)) then
            if (best_s == 0.or.dots(s) < dots(best_s)) best_s = s
          end if
        end do
        donor_to_site(i) = best_s
        used_site(best_s) = .true.
      end block
    end do

    do f = 1,nfrag
      unit_label(f) = fragments(f)%label
      unit_denticity(f) = fragments(f)%denticity
      unit_asym(f) = fragments(f)%asym
    end do

    call distinct_isomers(trim(geomname),nfrag,unit_label(1:nfrag),unit_denticity(1:nfrag), &
      & unit_asym(1:nfrag),max_bite_deg,isomers_out,niso,ok)
    if (.not.ok) then
      allocate (structures(3,nat,0))
      return
    end if

    call find_rotation_group(geom_sites,nsites,group,ngroup)

    !> build the input structure's own raw assignment, in the SAME
    !> (label,partner,end) representation distinct_isomers uses, so its
    !> canonical form can be matched against isomers_out and skipped
    input_raw(:,1:nsites) = 0
    do i = 1,ndonor_total
      f = donor_frag(i)
      s = donor_to_site(i)
      input_raw(1,s) = fragments(f)%label
      if (fragments(f)%denticity == 2) then
        !> find this donor's partner site: the OTHER donor of the same fragment
        block
          integer :: ii,other_donor_idx,psite
          other_donor_idx = 0
          do ii = 1,ndonor_total
            if (donor_frag(ii) == f.and.ii /= i) other_donor_idx = ii
          end do
          psite = donor_to_site(other_donor_idx)
          input_raw(2,s) = psite
          if (donor_atom(i) == fragments(f)%donors(1)) then
            input_raw(3,s) = 1
          else
            input_raw(3,s) = 2
          end if
        end block
      else
        input_raw(2,s) = 0
        input_raw(3,s) = 1
      end if
    end do
    call canonical_form(input_raw(:,1:nsites),nsites,group,ngroup,input_canon(:,1:nsites))

    allocate (structures(3,nat,niso))  !> at most niso-1 actually used; trimmed at the end
    nstruct = 0

    do k = 1,niso
      input_matches = all(isomers_out(:,:,k) == input_canon(:,1:nsites))
      if (input_matches) cycle

      new_coords = xyz_aa
      frag_of_label_consumed = 0

      !> walk sites in order, grouping partnered pairs, matching each
      !> site-group to the next not-yet-consumed same-label fragment (in
      !> original detection order) -- same-label fragments/site-groups
      !> are chemically interchangeable, so any such bijection yields
      !> the identical resulting molecule (see the Python reference's
      !> comment on this point)
      block
        logical :: site_done(nsites)
        integer :: lbl,partner,e_s,e_p,target_frag
        real(wp) :: old1(3),old2(3),new1(3),new2(3),newd(3)
        integer :: other_atoms(nat),nother,a
        real(wp) :: trial_coords(3,nat)
        real(wp),parameter :: trial_expand(4) = [0.0_wp,0.15_wp,0.3_wp,0.5_wp]
        integer :: itry,nclash,best_try
        real(wp) :: best_coords(3,nat)
        integer :: best_clash

        site_done = .false.
        do s = 1,nsites
          if (site_done(s)) cycle
          lbl = isomers_out(1,s,k)
          partner = isomers_out(2,s,k)

          !> find the next not-yet-consumed fragment with this label, in
          !> original detection order
          target_frag = 0
          do f = 1,nfrag
            if (fragments(f)%label == lbl.and.frag_of_label_consumed(f) == 0) then
              target_frag = f
              exit
            end if
          end do
          frag_of_label_consumed(target_frag) = 1
          site_done(s) = .true.

          nother = 0
          do a = 1,nat
            if (a == metal_idx) cycle
            if (any(fragments(target_frag)%atoms(1:fragments(target_frag)%natoms) == a)) cycle
            nother = nother+1
            other_atoms(nother) = a
          end do

          if (partner == 0) then
            !> denticity 1: free twist about the bond axis, keep the
            !> candidate with fewest clashes against everything else
            old1 = (xyz_aa(:,fragments(target_frag)%donors(1))-metal_pos)
            old1 = old1/norm2(old1)
            new1 = geom_sites(:,s)
            block
              real(wp) :: R0(3,3),Rtwist(3,3)
              integer :: tdeg,at_idx,ii
              real(wp) :: rel(3)
              R0 = align_vectors(old1,new1)
              best_clash = huge(1)
              do tdeg = 0,345,15
                Rtwist = matmul(rotmat_axis_angle(new1,tdeg*acos(-1.0_wp)/180.0_wp),R0)
                trial_coords = new_coords
                do ii = 1,fragments(target_frag)%natoms
                  at_idx = fragments(target_frag)%atoms(ii)
                  rel = xyz_aa(:,at_idx)-metal_pos
                  trial_coords(:,at_idx) = metal_pos+0.0_wp*new1+matmul(Rtwist,rel)
                end do
                nclash = clash_count(nat,at,trial_coords,fragments(target_frag)%natoms, &
                  & fragments(target_frag)%atoms(1:fragments(target_frag)%natoms),nother, &
                  & other_atoms(1:nother))
                if (nclash < best_clash) then
                  best_clash = nclash
                  best_coords = trial_coords
                end if
                if (best_clash == 0) exit
              end do
            end block
            new_coords = best_coords
          else if (partner > s) then
            !> denticity 2, process the pair once (when we hit the lower
            !> site index of the pair)
            site_done(partner) = .true.
            e_s = isomers_out(3,s,k); e_p = isomers_out(3,partner,k)
            if (e_s == 1) then
              new1 = geom_sites(:,s); new2 = geom_sites(:,partner)
            else
              new1 = geom_sites(:,partner); new2 = geom_sites(:,s)
            end if
            old1 = (xyz_aa(:,fragments(target_frag)%donors(1))-metal_pos)
            old1 = old1/norm2(old1)
            old2 = (xyz_aa(:,fragments(target_frag)%donors(2))-metal_pos)
            old2 = old2/norm2(old2)

            block
              real(wp) :: R(3,3),outward(3),rel(3)
              integer :: ii,at_idx
              R = exact_two_point_rotation(old1,new1,old2,new2)
              outward = new1+new2
              outward = outward/norm2(outward)
              best_clash = huge(1)
              do itry = 1,4
                trial_coords = new_coords
                do ii = 1,fragments(target_frag)%natoms
                  at_idx = fragments(target_frag)%atoms(ii)
                  rel = xyz_aa(:,at_idx)-metal_pos
                  trial_coords(:,at_idx) = metal_pos+trial_expand(itry)*outward+matmul(R,rel)
                end do
                nclash = clash_count(nat,at,trial_coords,fragments(target_frag)%natoms, &
                  & fragments(target_frag)%atoms(1:fragments(target_frag)%natoms),nother, &
                  & other_atoms(1:nother))
                if (nclash < best_clash) then
                  best_clash = nclash
                  best_coords = trial_coords
                end if
                if (best_clash == 0) exit
              end do
            end block
            new_coords = best_coords
          end if
        end do
      end block

      !> final whole-structure clash report
      block
        integer :: total_clash,a,other_atoms(nat),nother
        total_clash = 0
        do f = 1,nfrag
          nother = 0
          do a = 1,nat
            if (a == metal_idx) cycle
            if (any(fragments(f)%atoms(1:fragments(f)%natoms) == a)) cycle
            nother = nother+1
            other_atoms(nother) = a
          end do
          total_clash = total_clash+clash_count(nat,at,new_coords,fragments(f)%natoms, &
            & fragments(f)%atoms(1:fragments(f)%natoms),nother,other_atoms(1:nother))
        end do
        total_clash = total_clash/2
        if (total_clash > 0) then
          write (stdout,'(1x,a,i0,a,i0,a)') 'isomer_generate: WARNING ',total_clash, &
            & ' steric clash(es) in generated isomer ',k,' -- inspect before use'
        end if
      end block

      nstruct = nstruct+1
      structures(:,:,nstruct) = new_coords/autoaa   !> back to Bohr
    end do

    if (nstruct < niso) then
      block
        real(wp),allocatable :: tmp(:,:,:)
        allocate (tmp(3,nat,nstruct))
        tmp = structures(:,:,1:nstruct)
        deallocate (structures)
        allocate (structures(3,nat,nstruct))
        structures = tmp
      end block
    end if
    ok = .true.
  end subroutine generate_isomer_structures

end module isomer_generate_mod
