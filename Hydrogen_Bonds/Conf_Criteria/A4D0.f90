program criteria
!**********************************
!
! This fortran code is to read from file .fss_angs, .hbcase,
! and distinguish if the geometrical configuration around O*H is a perfect A4D0/A3D1.
! Modified from Hydrogen_Bonds/Length_Distribution/Codes/delta.f90
!
! Author: Lixin Zheng
! March 2015
!
implicit none
!
real, parameter         :: celldm = 23.5170     !dimension of the supercell 
real, parameter         :: convertBA= 0.52918     !convert Bohr into Angestrom
real, parameter         :: pi=4.D0*DATAN(1.D0)
REAL, PARAMETER         :: r_cov=1.24
character(len=40)       :: filename
character(len=20)       :: order
integer                 :: i,j,k,p,q,   &       !General index of atom,
                           iiO,iiH,     &       !the O* index
                           nsp(2),      &       !number of atoms in O and H
                           step,        &
                           stepstart, stepstop
integer                 :: step_t=0,iiO_t      !the step and oxygen index number from .trans file
integer                 :: total(4)
integer                 :: cs
integer                 :: state                !state=2: >0.5
                                                !state=4: <0.1
                                                
integer                 :: numION,numH
integer                 :: ncount=0
integer                 :: ierror
integer                 :: iO(15),iH(30),hnum(4),hnum_t(4)
real                    :: time
real                    :: rO(3,15),rH(3,30),rIO(3),rIH(3)
real                    :: rIO_t(3),rIH_t(3)
real                    :: bond(2)              !the calculated bond length
real                    :: sum(2)               !for the bond length when not transiting ions
real                    :: per                  !percentage of stable ions among all steps.
real                    :: sigma
real                    :: d(4)
! Local Variables
integer :: hb(4),iH_cov(4),hba,hbd
integer :: dummy,ctn,flag
integer :: pp,p1,p2
integer :: perfect_A4D0=0
real    :: angle(3),da(3)
real    :: rO_hb(3,4),rH_hb(3,4)
integer :: index_hb(4,2)
integer :: right_count=0
real    :: right_angle=pi/2
real    :: angle_thr=pi/18    ! 10 degree allowance of fluctration
!real    :: angle_thr=pi/18   ! 10 degree allowance of fluctration
real    :: eigen_angle

!
namelist /input/ cs, filename, stepstart, stepstop
!
!initialization
call init
!
!Open files
call files1
!
call main
!
call files2
!
!
! 
contains
!******************************************************************************************
!
!******************************************************************************************
subroutine init
!
implicit none
!read Namelist input for stdin
read(*,input)
!
!print Intro
write(*,*) '------------------------------------------'
write(*,*) '|          Find the Bond Length           |'
write(*,*) '------------------------------------------'
write(*,*) '1st solvation file : ', trim(filename)//'.fss_angs'
write(*,*) 'HB  file           : ', trim(filename)//'.hbcase'
write(*,*) 'A4D0 file          : ', trim(filename)//'.A4D0'
!
if (cs .eq. 2) numH=1
if (cs .eq. 3) numH=3
!
sum=0
total=0
!
end subroutine init
!******************************************************************************************
!
!******************************************************************************************
subroutine files1
implicit none
open(unit=3, file=(trim(filename)//'.fss_angs'), status='old')
open(unit=4, file=(trim(filename)//'.hbcase'), status='old')
open(unit=5, file=(trim(filename)//'.A4D0'), status='unknown')
end subroutine files1
!
subroutine files2
implicit none
close(3)
close(4)
close(5)
end subroutine files2 
!
!******************************************************************************************
subroutine main
!
implicit none
!
!************************************
!Jump out to stepstart
!if (stepstart .gt. 1) then
!  do k=1, stepstart-1
!    read(3,*) step,time,numION,nsp(1),nsp(2)
!    do i=1,nsp(1)+nsp(2)+1
!      read(3,*)
!    enddo
!    read(4,*) 
!  enddo
!endif
!************************************
!main loop
main_loop: do
  !
  flag=0
  right_count=0
  !if (ncount .gt. 10000) exit
  !
  read(3,*,iostat=ierror) step,time,numION,nsp(1),nsp(2)
  if (ierror .lt. 0) then
    write(*,*) '  End of File Reached'
    write(*,fmt='(1X, "  Total number of Samples: ", I7)' )(ncount)
    exit
  endif
  !**********************************
  read(4,*,iostat=ierror) step_t,time,dummy,dummy,hba,hbd,hb(1:4),dummy
  if (ierror .lt. 0) then
    write(*,*) '  End of File Reached'
    write(*,fmt='(1X, "  Total number of Samples: ", I7)' )(ncount)
    exit
  endif
  !**********************************
  !
  if (hba .eq. 4 .and. hbd .eq. 0) flag=1
  !if(flag .eq. 1) write(*,*) hba,hb(1:hba)
  !
  !**********************************
  !
  if (nsp(1) .eq. 0) then
    write(*,*) 'Error! nsp1=0.'
    exit 
  endif
  !
  !**********************************
  !Read in .fss file
  !**********************************
  q=0
  read(3,*) iiO,iiH
  !write(*,*) "rO_hb:"
  do i=1,nsp(1)
    !Read in oxygen positions
    read(3,*)  iO(i),rO(1:3,i),iH_cov(1:4)
    if (iO(i) .eq. iiO) then
      do j=1,3
        rIO(j)=rO(j,i)
      enddo 
      !if(flag .eq. 1) write(*,*) "rIO:", rIO(1:3)
    endif
    !======
    ! Find the 4 Hydrogen-bonded oxygen.
    if (flag .eq. 1) then
      ctn=0
      do p=1,4
        if (iH_cov(p) .eq. 0) exit
        if (iH_cov(p) .eq. hb(q+1)) then
          ! This O is HB-bonded to the ion
          ctn=1
          q=q+1
          exit
        endif
      enddo
      ! If this O is HB-bonded to the ion
      if (ctn .eq. 1) then
        do j=1,3
          rO_hb(j,q)=rO(j,i)
        enddo
        index_hb(q,1)=iO(i)
        index_hb(q,2)=hb(q)
        !if(flag .eq. 1) write(*,*) q, rO_hb(1:3,q)
      endif
    endif
  enddo !i=1,nsp1
  !
  if (q .gt. hba) write(*,*) "Error! Found more HB-bonded oxygen than indicated from *.hbcase!"
  !
  q=0
  !write(*,*) "rH_hb:"
  do i=1,nsp(2)
    !Read in hydrogen positions
    read(3,*)  iH(i),rH(1:3,i)
    if (iH(i) .eq. iiH) then
      do j=1,3
        rIH(j)=rH(j,i)
      enddo
    endif
    !======
    ! Find the 4 Hydrogen-bonded hydrogen.
    if (flag .eq. 1) then
      ctn=0
      pp=0
      do p=1,4
        if (iH(i) .eq. index_hb(p,2)) then
          ctn=1
          pp=p
        endif
      enddo
      if (ctn .eq. 1) then
        ! the order of these 4 donated Hydrogen is according to original 1-127 index.
        do j=1,3
          rH_hb(j,pp)=rH(j,i)
        enddo
      endif
    endif
  enddo !i=1,nsp2
  !**********************************
  !
  if (flag .eq. 0) cycle main_loop
  !
  !**********************************
  !Core Contents
  do q=1,4
    pp=1
    do p1=1,3
      if (p1 .eq. q) cycle
      do p2=p1+1,4
        if (p2 .eq. q) cycle
        call get_angle(rH_hb(1:3,q),rH_hb(1:3,p1),rH_hb(1:3,p2),angle(pp))
        !write(*,*) "q,p1,p2:",q,p1,p2,angle(pp)
        da(pp)=abs(angle(q)-pi/2)
        pp=pp+1
      enddo
    enddo
    if (pp .ne. 4) then
      if(flag .eq. 1) write(*,*) "error! pp .ne. 3!"
      exit main_loop
    endif
    ! ?
    pp=1
    !if(da(2) .gt. da(1)) then
    !  pp=2
    !  da(1)=da(2)
    !endif
    !if(da(3) .gt. da(1)) then
    !  pp=3
    !  da(1)=da(2)
    !endif
    if (angle(2) .gt. angle(1)) then
      pp=2
      angle(1)=angle(2)
    endif
    if (angle(3) .gt. angle(1)) then
      pp=3
      angle(1)=angle(3)
    endif
    !
    if(q .eq. 1) then
      if (pp .eq. 1) order=" 1  2  4  3"
      if (pp .eq. 2) order=" 1  2  3  4"
      if (pp .eq. 3) order=" 1  3  2  4"
      !write(*,*) order
    endif
    !
    !write(*,*) "The correct angle of hb",q,"is",pp
    !write(*,*) q,angle(1)
    eigen_angle=angle(1)
    if (abs(eigen_angle-right_angle) .lt. angle_thr) right_count=right_count+1
  enddo !q
  if (right_count .eq. 4) then
    perfect_A4D0=perfect_A4D0+1
    write(5,fmt='(I10,F16.8,8I6,3x,a20)') step,time, index_hb(1:4,1), index_hb(1:4,2), order
    !write(5,*) "O",rIO(1:3)
    !do p=1,4
    !  write(5,*) "H",rH_hb(1:3,p)
    !enddo
    !write(*,*) "Find perfect A4D0 at",step, time
  endif
  !
  !Core Contents ends
  !**********************************
  !
  ncount=ncount+1

enddo main_loop

write(*,*) "Percentage of perfect A4D0:"
write(*,*) perfect_A4D0,real(perfect_A4D0)/real(ncount)*100

end subroutine main
!
!******************************************************
!******************************************************
subroutine get_distance(pos1,pos2,dist)

implicit none
!
real, intent(in)    :: pos1(3),pos2(3)
real, intent(inout) :: dist
real                :: delta(3)
real, parameter     :: pcell=23.5170*convertBA
integer             :: i
!
do i=1,3
  delta(i)=pos1(i)-pos2(i)
  delta(i)=delta(i)-nint(delta(i)/pcell)*pcell
enddo
dist=sqrt(delta(1)**2+delta(2)**2+delta(3)**2)
!
return
end subroutine get_distance
!******************************************************
subroutine get_angle(pos1,pos12,pos13,a)

implicit none
!
real, intent(in)    :: pos1(3),pos12(3),pos13(3)
real, intent(inout) :: a
real                :: delta(3,2),dist1,dist2,dot
real, parameter     :: pcell=23.5170*convertBA
integer             :: i
!
do i=1,3
  delta(i,1)=pos1(i)-pos12(i)
  delta(i,1)=delta(i,1)-nint(delta(i,1)/pcell)*pcell
  !
  delta(i,2)=pos1(i)-pos13(i)
  delta(i,2)=delta(i,2)-nint(delta(i,2)/pcell)*pcell
enddo

dot=delta(1,1)*delta(1,2)+&
    delta(2,1)*delta(2,2)+&
    delta(3,1)*delta(3,2)
dist1=sqrt(delta(1,1)**2+delta(2,1)**2+delta(3,1)**2)
dist2=sqrt(delta(1,2)**2+delta(2,2)**2+delta(3,2)**2)

a=acos(dot/(dist1*dist2))
!
return
end subroutine get_angle
!******************************************************************************************
end program criteria
