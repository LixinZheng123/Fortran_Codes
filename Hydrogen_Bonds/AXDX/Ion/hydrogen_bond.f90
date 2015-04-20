program hydrogen_bond
!**********************************
!
!This fortran file is to read from file .ion, .fss, 
!and calculate the total AXDX percentage from ion to first solvation shell atoms.
!Author: Lixin Zheng
!Sep 2014
!
implicit none
!
INTEGER, PARAMETER      :: DP = selected_real_kind(14,200) !double-precision kind
real(DP), parameter     :: pi=4.d0*atan(1.d0)
real(DP), parameter     :: celldm = 23.5170     !dimension of the supercell 
real(DP), parameter     :: convertBA= 0.52918     !convert Bohr into Angestrom
character(len=40)       :: filename
integer                 :: i,j,k,p,q,      &    !General index of atom,
                           iO(15),iH(30),  &    !the index of first solvation shell (FSS)
                           step, step_t=0, step_p,   &
                           stepstop
integer                 :: nsp(2)
integer                 :: total_accept=0,total_donate=0 
integer                 :: cs
integer                 :: state                
integer                 :: readstate                
integer                 :: countS(4)
integer                 :: numION,numH
integer                 :: ncount=0
integer                 :: ierror
integer                 :: iiO, iiH(3), iiO_t, iiH_t
integer                 :: num_donate,num_accept,HBnum(2),HBIndex(6,2)
integer                 :: AXDX(5,5,4)                    !Accept, donate, state
integer                 :: cov_h(4,15)
integer                 :: hbcase
real(DP)                :: time
real(DP)                :: rO(3,15), rH(3,30)
real(DP)                :: rIO(3), rIO_t(3), rIH(3,3), rIH_t(3)
real(DP)                :: delta(3,2)
real(DP)                :: d(2)
real(DP)                :: theta
real(DP)                :: rHB
real(DP)                :: percent(5,5)
!
namelist /input/ cs, filename, rHB, stepstop
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
write(*,*) '1st Shell file       : ', trim(filename)//'.fss_angs'
write(*,*) 'Output file          : ', trim(filename)//'.hbcase'
!
if (cs .eq. 2) numH=1
if (cs .eq. 3) numH=3
!
AXDX=0
countS=0
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
    open(unit=3, file=(trim(filename)//'.fss_angs'), status='old')
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
  !Read .fss
  read(3,*,iostat=ierror) step,time,numION,nsp(1),nsp(2)
  if (ierror .lt. 0) then
    write(*,*) '  End of File Reached'
    write(*,fmt='(1X, "  Total number of Samples: ", I7)' )(ncount)
    exit
  endif
  !**********************************
  !Error message
  if (nsp(1) .eq. 0) then
    write(*,*) 'Error! nsp1=0.'
    exit
  endif
  !**********************************
  !Read .fss
  read(3,*) iiO,iiH(1:numH)
  do i=1,nsp(1)
    read(3,*)  iO(i),rO(1:3,i),cov_h(1:4,i)
    if (iO(i) .eq. iiO) then
      do j=1,3
        rIO(j)=rO(j,i)
      enddo
    else if (iO(i) .eq. iiO_t) then
      do j=1,3
        rIO_t(j)=rO(j,i)
      enddo
    endif
  enddo
  !
  q=0
  HBIndex=0
  !
  do i=1,nsp(2)
    read(3,*)  iH(i),rH(1:3,i)
    if (iH(i) .eq. iiH(q+1) .and. q .lt. numH) then
      q=q+1
      do j=1,3
        rIH(j,q)=rH(j,i)
      enddo
    endif
    if (iH(i) .eq. iiH_t) then
      do j=1,3
        rIH_t(j)=rH(j,i)
      enddo
    endif
  enddo
  !********************************************
  !
  num_accept=0
  num_donate=0
  !
  !===============count the ACCEPTING Hydrogen Bonds=================
  do p=1,nsp(1)
    call get_distance(rO(1:3,p),rIO(1:3),d(1))
    if (d(1) .lt. rHB) then
      do q=1,4
        if (cov_h(q,p) .eq. 0) exit
        do i=1, nsp(2)
          if (cov_h(q,p) .eq. iH(i)) then
            call get_distance(rH(1:3,i),rO(1:3,p),d(2))
            call get_theta(rO(1:3,p),rIO(1:3),rH(1:3,i),d(1:2),theta)
            if (theta .le. pi/6) then 
              num_accept=num_accept+1
              if (num_accept .gt. 5) then
                write(*,'("Error! HB of O* equals to ",I1,"at step= ", I8)') &
                     num_accept, step
              endif
              HBIndex(num_accept,1)=iH(i)
              !write(*,*) "Accepting from iH:", i
            endif
          endif
        enddo
      enddo
    endif
  enddo
  !===============count the DONATING Hydrogen Bonds=================
  do q=1,numH
    do p=1,nsp(1)
      call get_distance(rO(1:3,p),rIO(1:3),d(1))
      if (d(1) .lt. rHB) then
        call get_distance(rIH(1:3,q),rIO(1:3),d(2))
        call get_theta(rIO(1:3),rO(1:3,p),rIH(1:3,q),d(1:2),theta)
        if (theta .le. pi/6) then
          num_donate=num_donate+1
          HBIndex(num_donate,2)=iO(p)
          !write(*,*) "Donating from iO to iiH:", i, iiH(q)
        endif
      endif
    enddo
  enddo
  !================================================================
  total_donate=total_donate+num_donate
  total_accept=total_accept+num_accept
  HBnum(1)=num_accept
  HBnum(2)=num_donate
  !write(*,*) "At ncount=", ncount+1
  !write(*,*) "number of donating/accepting:      ", num_donate, num_accept
  !write(*,*) "number of total donating/accepting:", total_donate, total_accept
  !================================================================
  !For OH-,  the important ones are A4D0*, A3D0*, A3D1*, A2D2
  !For H3O+, the important ones are A0D3*, A1D3, A1D2, A2D2
  AXDX(num_accept+1,num_donate+1,state)=AXDX(num_accept+1,num_donate+1,state)+1
  !
  ncount=ncount+1
  countS(state)=countS(state)+1
  !
  !Modified 20140902
  hbcase=0
  if (num_accept .eq. 3 .and. num_donate .eq. 0) hbcase=1
  if (num_accept .eq. 3 .and. num_donate .eq. 1) hbcase=2
  if (num_accept .eq. 4 .and. num_donate .eq. 0) hbcase=3
  if (num_accept .eq. 4 .and. num_donate .eq. 1) hbcase=4
  !write(7,fmt='(I10,3X,3I5,3X,2I2,3X,4I5,3X,4I5)') &
  write(7,fmt='(I10,3X,F10.3,3X,2I5,3X,6I5,3X,I5)') &
      step, time, HBnum(1:2),HBIndex(1:6,1),HBIndex(1,2)
  !Modified 20140902 ends
  !
enddo
!
!
!
!********************writing out*******************
!
write(*,*)
write(*,*)
write(*,fmt='(5X,"These are the cases when percentage is greater than 1%:")')
do i=1,5
  do j=1,5
    percent(i,j)=(AXDX(i,j,1)+AXDX(i,j,2)+AXDX(i,j,3)+AXDX(i,j,4))/real(ncount)*100
    if (percent(i,j) .gt. 1) write(*,fmt='(5X,"A",I1,"D",I1,2X,F5.1)') i-1,j-1,percent(i,j)
  enddo
enddo
write(*,*)
write(*,fmt='(5X,"HB accepting/donating percentage in the concerning cases:")')
write(*,*)
if (cs .eq. 2) then
  write(*,fmt='(5X,"A3D0",4X,"A3D1",4X,"A4D0",4X,"A4D1")')
  write(*,fmt='(4X,F5.1,3X,F5.1,3X,F5.1,3X,F5.1)') percent(4,1), &
                                                   percent(4,2), &
                                                   percent(5,1), &
                                                   percent(5,2)
  !write(9,fmt='(13X,"A3D0",4X,"A3D1",4X,"A4D0",4X,"A4D1")')
  !write(9,fmt='("Sum",9X,F5.1,3X,F5.1,3X,F5.1,3X,F5.1)') percent(4,1), &
  !                                                       percent(4,2), &
  !                                                       percent(5,1), &
  !                                                       percent(5,2)
  !do state=1,4
  !  write(9,fmt='(I1,1X,"(",F4.1,"%)",3X,F5.1,3X,F5.1,3X,F5.1,3X,F5.1)') state, countS(state)/real(ncount)*100, &
  !                                                                    AXDX(4,1,state)/real(countS(state))*100, &
  !                                                                    AXDX(4,2,state)/real(countS(state))*100, &
  !                                                                    AXDX(5,1,state)/real(countS(state))*100, &
  !                                                                    AXDX(5,2,state)/real(countS(state))*100
  !enddo
!
else if (cs .eq. 3) then
  write(*,fmt='(5X,"A0D3",4X,"A1D3",4X,"A2D3")')
  write(*,fmt='(4X,F5.1,3X,F5.1,3X,F5.1)') percent(1,4), &
                                           percent(2,4), &
                                           percent(3,4)
  !write(9,fmt='(13X,"A0D3",4X,"A1D3",4X,"A2D3")')
  !write(9,fmt='("Sum",9X,F5.1,3X,F5.1,3X,F5.1)') percent(1,4), &
  !                                               percent(2,4), &
  !                                               percent(3,4)
  !do state=1,4
  !  write(9,fmt='(I1,1X,"(",F4.1,"%)",3X,F5.1,3X,F5.1,3X,F5.1)') state, countS(state)/real(ncount)*100, &
  !                                                            AXDX(1,4,state)/real(countS(state))*100, &
  !                                                            AXDX(2,4,state)/real(countS(state))*100, &
  !                                                            AXDX(3,4,state)/real(countS(state))*100
  !enddo
endif
!
!write(9,*)
write(*,*)
write(*,fmt='(5X,"Average accepting:  ",F4.2)') total_accept/real(ncount)
!write(9,fmt='("Average accepting:  ",F4.2)') total_accept/real(ncount)
write(*,fmt='(5X,"Average donating:   ",F4.2)') total_donate/real(ncount)
!write(9,fmt='("Average donating:   ",F4.2)') total_donate/real(ncount)
write(*,*)
write(*,*)
!!
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
