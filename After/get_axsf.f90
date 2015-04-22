program get_axsf
!
!This fortran file is to find first solvation shell (FSS) of a defect from a MD simulation trajectory.
!It is used for visualization and selection using VMD.
!Author: Lixin Zheng
!May 2014
!
implicit none
!
real, parameter         :: celldm = 23.5170     !dimension of the supercell 
real, parameter         :: convertBA= 0.529     !Transfer Bohr to Angestrom9
REAL, PARAMETER         :: r_FSS=3.5
REAL, PARAMETER         :: r_cov=1.24

character(len=40)       :: filename
character               :: dummyC
integer                 :: i,iO,iH,     &       !
                           j,k,p,q,q1,  &       !General index of atom,
                           iiO(5),      &       !the O* index
                           iiO_record(5),&      !the various O* index to be recorded in .fss
                           iFSS_O(15),iFSS_H(30),&!the index of first solvation shell
                           nsp1, nsp2,  &       !number of atoms in O and H
                           new_nsp1, new_nsp2,&
                           stepstop, readstep, readstepp, e_step
integer                 :: readstep_p=0
integer                 :: cs, numH
integer                 :: start=0
integer                 :: numION               !number of Ion
integer                 :: numW                 !number of water
integer                 :: iiH(5,4)             !(numH:numION)
integer                 :: iiO_p(2),iiHp(4)
integer                 :: iiO_p_record(2)      !The iiO that has been recorded to *.trans
integer                 :: iiH_real
integer                 :: error
real                    :: dr(5)
real                    :: rO(3,64),rH(3,128)
integer                 :: ierror,ncount=0
integer                 :: h_num(4,15)           !For each Oxygen in FSS
                                                !there should be 1 or 2 hydrogen linked to it
real                    :: d,           &       !distance
                           tau_O(3,15), tau_H(3,30)   !position of each FSS atoms in angstrom
real                    :: time
real                    :: r_iiO(3,5)
real                    :: per
integer                 :: i_anim=0
integer                 :: stopp
integer                 :: dstep
!
namelist /input/ filename, start, stopp, dstep, nsp1, nsp2
!
!initialization
call init
!
!Open files
call files1
!
!main loop
call main
!
call files2
!
!
! 
contains
!
!=========================================================================================
!=========================================================================================
!
subroutine init
!
implicit none
!read Namelist input for stdin
read(*,input)
!
!print Intro
write(*,*) '--------------------------------------------------'
write(*,*) '|     Find the Ion and First Solvation Shell      |'
write(*,*) '--------------------------------------------------'
write(*,*) '1st solvation file : ', trim(filename)//'.fss_angs'
write(*,*) 'Output file        : ', trim(filename)//'.axsf'

!
end subroutine init
!
!=========================================================================================
!=========================================================================================
!
subroutine files1
implicit none
open(unit=11, file=(trim(filename)//'.fss_angs'), status='unknown')
open(unit=13, file=(trim(filename)//'.axsf'), status='unknown')
end subroutine files1
!
subroutine files2
implicit none
close(11)
close(13)
end subroutine files2 
!
!=========================================================================================
!=========================================================================================
!
subroutine main
!
implicit none
!
!************************************
write(13,*) "ANIMSTEPS ",((stopp-start)/dstep+1)
write(13,*) "CRYSTAL"
write(13,*) "PRIMVEC"
write(13,*) "12.444660414    0.000000000    0.000000000"
write(13,*) " 0.000000000   12.444660414    0.000000000"
write(13,*) " 0.000000000    0.000000000   12.444660414"
!main loop
do 
  !
  read(11,*) readstep, time, numION, new_nsp1,new_nsp2
  if (readstep .gt. stopp) exit
  if (readstep .lt. start) then
    do q=1,new_nsp1+new_nsp2+1
      read(11,*)
    enddo
  else !if (readstep .ge. start) then
    !
    i_anim=i_anim+1
    if (i_anim .gt. ((stopp-start)/dstep+1)) write(*,*) "Error!"
    write(13,*) "PRIMCOORD   ", i_anim
    write(13,*) nsp1+nsp2, "      1"
    read(11,*) ! iiO(1),iiH(1:numH,1)
    !
    if (new_nsp1 .eq. nsp1) then
    ! The normal case
      do p=1,new_nsp1
        read(11,*) iFSS_O(p),tau_O(1:3,p)
        write(13,*) "O   ",tau_O(1:3,p)
      enddo
      do q=1,new_nsp2
        read(11,*) iFSS_H(q),tau_H(1:3,q)
        write(13,*) "H   ",tau_H(1:3,q)
      enddo
      !
    else ! new_nsp1 .ne. nsp1
      write(*,*) "Error! Please examine the nsp1 input!"
    endif
  endif
enddo

!
end subroutine main
end program get_axsf
