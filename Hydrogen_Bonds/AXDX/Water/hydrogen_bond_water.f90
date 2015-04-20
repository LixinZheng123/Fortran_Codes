program hydrogen_bond
!**********************************
!
! This fortran file is to read from file .pos, 
! and calculate the HB situation of every H2O,
! Also, it will print out the AXDX percentage.
! Modified from Analysis_Tools/Hydrogen_Bonds/AXDX/Ion/Codes/hydrogen_bond.f90
! Author: Lixin Zheng
! Mar 2015
!
implicit none
!
INTEGER, PARAMETER      :: DP = selected_real_kind(14,200) !double-precision kind
real(DP), parameter     :: pi=4.d0*atan(1.d0)
real,parameter          :: r_cov=1.24
real,parameter          :: rHB=3.5
real(DP), parameter     :: celldm = 23.5170     !dimension of the supercell 
real(DP), parameter     :: convertBA= 0.52918     !convert Bohr into Angestrom
character(len=40)       :: filename
integer                 :: i,i1,i2,j,k,p,q,      &    !General index of atom,
                           iO(15),iH(30),  &    !the index of first solvation shell (FSS)
                           step, step_t=0, step_p,   &
                           stepstart,stepstop
integer                 :: nsp(2)
integer                 :: total_accept=0,total_donate=0 
integer                 :: cs
!integer                 :: state                
integer                 :: readstate                
integer                 :: countS(4)
integer                 :: numION,numH
integer                 :: ncount=0
integer                 :: ierror
integer                 :: iiO, iiH(3) !, iiO_t, iiH_t
integer                 :: num_donate,num_accept
integer                 :: AXDX(5,5,4)                    !Accept, donate, state
integer                 :: cov_h(3,70)                    ! Covalent-bonded hydrogens
integer                 :: hbcase
integer                 :: HBIndex(5,2,64)
integer                 :: HBnum(2,64)
real(DP)                :: time
real(DP)                :: sigma
real(DP)                :: rO(3,70), rH(3,140)
real(DP)                :: rIO(3), rIO_t(3), rIH(3,3), rIH_t(3)
real(DP)                :: delta(3,2)
real(DP)                :: d(2)
real(DP)                :: theta
real(DP)                :: percent(5,5)
real(DP)                :: pct_11, pct_12, pct_21, pct_22, pct_32
!
namelist /input/ cs, filename, stepstart,stepstop
!
!initialization
call init
!
!Open files
call files(1)
!
call main
!
call files(2)
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
write(*,*) 'Input file       : ', trim(filename)//'.pos '
write(*,*) 'Output file      : ', trim(filename)//'.hbcase'
!write(*,*) 'Time-correlation file: ', trim(filename)//'.'
!
nsp(1)=64
nsp(2)=128
!
!
end subroutine init
!******************************************************************************************
!
!******************************************************************************************
subroutine files(io)
!
implicit none
!
integer,intent(in) :: io
!
select case(io)
  case(1)
    open(unit=3, file=(trim(filename)//'.pos'), status='old')
    open(unit=7, file=(trim(filename)//'.hbcase'), status='unknown')
  case(2)
    close(3)
    close(7)
  case default
end select

end subroutine files
!
!******************************************************************************************
!
!******************************************************************************************
subroutine main
!
implicit none
!
!************************************
!main loop
do
  !********************************************
  !if (ncount .gt. 100) exit
  !Read .pos
  read(3,*,iostat=ierror) step,time
  if (ierror .lt. 0) then
    write(*,*) '  End of File Reached'
    write(*,fmt='(1X, "  Total number of Samples: ", I7)' )(ncount)
    exit
  endif
  write(7,*) step,time
  !**********************************
  !Read .pos
  do i=1,nsp(1)
    read(3,*) rO(1:3,i)
    ! Change the unit to Angstrom
    do j=1,3
      rO(j,i)=rO(j,i)*convertBA
    enddo
  enddo
  do i=1,nsp(2)
    read(3,*) rH(1:3,i)
    ! Change the unit to Angstrom
    do j=1,3
      rH(j,i)=rH(j,i)*convertBA
    enddo
  enddo
  !********************************************
  !
  num_accept=0
  num_donate=0
  cov_h=0
  HBIndex=0
  !
  !Record the covalent bond between O and H
  do i1=1,nsp(1)
    q=0
    do i2=1,nsp(2)
      call get_distance(rO(1:3,i1),rH(1:3,i2),d(1))
      if (d(1) .lt. r_cov) then
        q=q+1
        cov_h(q,i1)=i2
      endif
    enddo
    if (q .ne. 2) then
      write(*,'("Error! Number of covalen-bonded H around O ",I2," &
           equals to ",I1," at step= ", I8)') i1,q,step
      exit
    endif
  enddo
  !===============count the ACCEPTING Hydrogen Bonds=================
  do i=1,nsp(1)
    p=0
    do i1=1,nsp(1)
      if (i .eq. i1) cycle
      call get_distance(rO(1:3,i),rO(1:3,i1),d(1))
      if (d(1) .gt. rHB) cycle
      !write(*,*) i1, "is inside the fss of",i
      !
      ! Should we include these with non-standard water structure?
      ! If the covalent bond of O(i1) is more than 2, we'll skip this i1.
      ! If the covalent bond of O(i1) is less than 2, we'll skip this i1.
      if (cov_h(3,i1) .ne. 0 .or. cov_h(2,i1) .eq. 0) cycle
      !
      do q=1,2
        call get_distance(rO(1:3,i1),rH(1:3,cov_h(q,i1)),d(2))
        call get_theta(rO(1:3,i1),rO(1:3,i),rH(1:3,cov_h(q,i1)),d(1:2),theta)
        !write(*,*) d(2),theta
        if (theta .le. pi/6) then
          p=p+1
          if (p .gt. 4) then
            write(*,'("Error! HB of O ",I2,"equals to ",I1,"at step= ", I8)') &
                 i, p, step
          endif
          HBIndex(p,1,i)=cov_h(q,i1)
          !write(*,*) "O #",i,", HB accepted",HBIndex(p,1,i)
        endif ! if (theta .le. pi/6)
      enddo ! do q=1,2
      !
    enddo ! do i1=1,nsp(1)
    !
    !num_accept=p
    !write(9,*) i,cov_h(1:2,i),hb_accept(1:4)
    HBnum(1,i)=p
  enddo ! do i=1,nsp(1)
  !
  !===============count the DONATING Hydrogen Bonds=================
  do i=1,nsp(1)
    p=0 
    do i1=1,nsp(1)
      if (i .eq. i1) cycle
      if (cov_h(p+1,i) .eq. 0) exit
      do q=1,HBnum(1,i1)
        if (cov_h(p+1,i) .eq. HBIndex(q,1,i1)) then
          p=p+1
          HBIndex(p,2,i)=i1
        endif
      enddo ! do q=1,4
    enddo ! do i1=1,nsp(1)
    HBnum(2,i)=p
    write(7,fmt='(I2,3X,3I5,3X,2I2,3X,4I5,3X,4I5)') &
         i,cov_h(1:3,i),HBnum(1:2,i),HBIndex(1:4,1,i),HBIndex(1:4,2,i)
  enddo ! do i=1,nsp(1)
  !num_donate
  ncount=ncount+1      
enddo

end subroutine main
!
!******************************************************
!Apply periodic boundary conditions to wrap 
! coordinates around the origin
!******************************************************
subroutine get_distance(pos1,pos2,dist)

implicit none
!
real(DP), intent(in)    :: pos1(3),pos2(3)
real(DP), intent(inout) :: dist
real(DP)                :: delta(3)
real(DP), parameter     :: pcell=23.5170*convertBA
integer                 :: i
!
do i=1,3
  delta(i)=pos1(i)-pos2(i)
  delta(i)=delta(i)-nint(delta(i)/pcell)*pcell
enddo
dist=sqrt(delta(1)**2+delta(2)**2+delta(3)**2)
!
return
end subroutine get_distance

subroutine get_theta(pos1,pos2,pos3,pdist,angle)
!Try to find the angle of line 1-2, and lind 1-3.
implicit none

real(DP), intent(in)    :: pos1(3),pos2(3),pos3(3)
real(DP), intent(in)    :: pdist(2)  !pdist(1) is r(1-2),pdist(2) is r(1-3)
real(DP), intent(inout) :: angle
real(DP)                :: delta(3,2)
real(DP), parameter     :: pcell=23.5170*convertBA
integer                 :: i
!
do i=1,3
  delta(i,1)=pos1(i)-pos2(i)
  delta(i,1)=delta(i,1)-nint(delta(i,1)/pcell)*pcell
  delta(i,2)=pos1(i)-pos3(i)
  delta(i,2)=delta(i,2)-nint(delta(i,2)/pcell)*pcell
enddo
angle=acos((delta(1,1)*delta(1,2)+delta(2,1)*delta(2,2)+delta(3,1)*delta(3,2))/(pdist(1)*pdist(2)))

return
end subroutine get_theta
!******************************************************************************************
end program hydrogen_bond
