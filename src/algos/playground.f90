!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2022 Philipp Pracht
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
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!> Implementation of whatever, for testing implementations
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!========================================================================================!
!> Input/Output:
!>  env  -  crest's systemdata object
!>  tim  -  timer object
!>-----------------------------------------------
subroutine crest_playground(env,tim)
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd
  use isomer_enumerate_mod
  use isomer_generate_mod
  implicit none
  logical,external :: isTMetal
  integer,allocatable :: iso_out(:,:,:)
  integer :: iso_n
  logical :: iso_ok
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich
  logical :: pr,wr
!========================================================================================!
  type(calcdata) :: calc
  real(wp) :: accuracy,etemp

  integer :: V,maxgen
  integer,allocatable :: A(:,:)
  logical,allocatable :: rings(:,:)
  integer,allocatable :: tmp(:)
  logical :: connected,fail,doreturn

  real(wp) :: energy
  real(wp),allocatable :: grad(:,:),geo(:,:),csv(:,:),q(:)

!========================================================================================!
  call tim%start(14,'Test implementation')
!========================================================================================!
  !call system('figlet welcome')
  write (*,*) "              _                          "
  write (*,*) "__      _____| | ___ ___  _ __ ___   ___ "
  write (*,*) "\ \ /\ / / _ \ |/ __/ _ \| '_ ` _ \ / _ \"
  write (*,*) " \ V  V /  __/ | (_| (_) | | | | | |  __/"
  write (*,*) "  \_/\_/ \___|_|\___\___/|_| |_| |_|\___|"
  write (*,*)

!>--- TEMPORARY: distinct_isomers self-checks (mirrors tools/isomer_enumerate.py --self-check)
  call distinct_isomers('octahedral',6,[1,1,1,2,2,2],[1,1,1,1,1,1], &
    & [.false.,.false.,.false.,.false.,.false.,.false.],100.0_wp,iso_out,iso_n,iso_ok)
  write (*,'(1x,A,I0,A)') 'fac/mer (octahedral 3A3B): found ',iso_n,' (expected 2)'
  deallocate (iso_out)

  call distinct_isomers('octahedral',3,[1,1,1],[2,2,2],[.false.,.false.,.false.],95.0_wp,iso_out,iso_n,iso_ok)
  write (*,'(1x,A,I0,A)') 'Bailar symmetric bidentate: found ',iso_n,' (expected 2)'
  deallocate (iso_out)

  call distinct_isomers('octahedral',3,[1,1,1],[2,2,2],[.true.,.true.,.true.],95.0_wp,iso_out,iso_n,iso_ok)
  write (*,'(1x,A,I0,A)') 'Bailar ASYMMETRIC bidentate: found ',iso_n,' (must be > symmetric''s 2)'
  deallocate (iso_out)

  call distinct_isomers('square_planar',4,[1,1,2,2],[1,1,1,1], &
    & [.false.,.false.,.false.,.false.],100.0_wp,iso_out,iso_n,iso_ok)
  write (*,'(1x,A,I0,A)') 'cis/trans (square_planar 2A2B): found ',iso_n,' (expected 2)'
  deallocate (iso_out)

  call distinct_isomers('tetrahedral',4,[1,1,2,2],[1,1,1,1], &
    & [.false.,.false.,.false.,.false.],100.0_wp,iso_out,iso_n,iso_ok)
  write (*,'(1x,A,I0,A)') 'tetrahedral 2A2B (should have NO cis/trans): found ',iso_n,' (expected 1)'
  deallocate (iso_out)

  call distinct_isomers('trigonal_bipyramidal',5,[1,1,2,2,2],[1,1,1,1,1], &
    & [.false.,.false.,.false.,.false.,.false.],100.0_wp,iso_out,iso_n,iso_ok)
  write (*,'(1x,A,I0,A)') 'TBP 2A3B: found ',iso_n,' (informational, python got 3)'
  deallocate (iso_out)
  write (*,*)

!>--- TEMPORARY: generate_isomer_structures end-to-end test on the loaded input
  block
    real(wp),allocatable :: structs(:,:,:)
    integer :: nstruct,ii
    character(len=24) :: geomname
    logical :: genok
    type(coord) :: outmol
    character(len=40) :: fname
    integer :: metal_idx,ii2
    call env%ref%to(mol)
    metal_idx = 1
    do ii2 = 1,mol%nat
      if (isTMetal(mol%at(ii2))) then
        metal_idx = ii2
        exit
      end if
    end do
    write (*,'(1x,A,I0)') 'auto-detected metal atom index: ',metal_idx
    call generate_isomer_structures(mol%nat,mol%at,mol%xyz,metal_idx,100.0_wp, &
      & structs,nstruct,geomname,genok)
    write (*,'(1x,A,L1,A,A,A,I0)') 'generate_isomer_structures: ok=',genok, &
      & ' geometry=',trim(geomname),' nstruct=',nstruct
    do ii = 1,nstruct
      outmol = mol
      outmol%xyz = structs(:,:,ii)
      write (fname,'(A,I0,A)') 'generated_isomer_',ii,'.xyz'
      call outmol%write(trim(fname))
      write (*,'(1x,A,A)') '  wrote ',trim(fname)
    end do
  end block
  write (*,*)
!========================================================================================!
  call env%ref%to(mol)
  write (*,*)
  write (*,*) 'Input structure:'
  call mol%append(stdout)
  write (*,*)
!!========================================================================================!
!
!  allocate (grad(3,mol%nat),source=0.0_wp)
!  call env2calc(env,calc,mol)
!  calc%calcs(1)%rdwbo = .true.
!  call calc%info(stdout)
!
!  call engrad(mol,calc,energy,grad,io)
!  call calculation_summary(calc,mol,energy,grad)
!========================================================================================!

  allocate(mol%gradient(3,mol%nat), source=1.0_wp)
  call mol%write('dummy.extxyz')


  call molnew%open("dummy.extxyz")
  call molnew%write("dummy2.extxyz")


  block
    type(coord),allocatable :: structures(:)
    integer :: nall
    call rdensemble(env%inputcoords,nall,structures)
    write(*,*) nall,'structures read from ',env%inputcoords
    call wrensemble('dummyensemble.xyz',nall,structures)
  end block
!========================================================================================!
  call tim%stop(14)
  return
end subroutine crest_playground
