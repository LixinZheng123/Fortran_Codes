program mean_square_displacement
!
implicit none
!
INTEGER, PARAMETER      :: DP = selected_real_kind(14,200) !double-precision kind
real(DP), parameter     :: celldm = 23.5170,&   !dimension of the supercell 
                           convertBA = 0.529    !Transfer Bohr to Angestrom
real(DP)                :: cutoff1=0.7*celldm   !
real(DP)                :: cutoff2=0.25*celldm  !
character(len=40)       :: filename 
integer                 :: cs
integer                 :: i,ip,j,p,    &       !General index of atom,
                           k,k0,dk,&            !steps
                           n,           &       !For counting numbers and normalize
                           nsp(2),      &       !number of atoms in O and H
                           numION,      &       !Number of ion in each timestep
                           dstep,       &       !step difference from *.ion
                           stepstart
integer                 :: dummyI
integer                 :: iO,iiO
integer                 :: ncount=0
integer                 :: ierror
integer                 :: last_dk              ! From 0 to last_dk is the range of msd output.
integer                 :: average_box(3)
integer                 :: e_step(5), e_iiO(5)
real(DP)                :: sum,         &
                           msd,         &       !mean square displacement
                           tot_sd,      &       !square displacement before nomalization
                           d2,          &
                           dummyR               !a dummy variable
real(DP)                :: tmp_dr                  
real(DP)                :: dr(3)
real(DP)                :: amass(2), pmass(2)
real(DP),allocatable    :: time(:),     &      
                           rIO(:,:), drp(:,:)
integer,allocatable     :: box(:,:),step(:)
!
namelist /input/ filename, stepstart, dstep, last_dk, cs
!
!initialization
call init
!
!Open files
call files1
!
!read data from *pos
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
!
!read Namelist input for stdin
read(*,input)
!
!print Intro
write(*,*) '------------------------------------------'
write(*,*) '|  Mean Square Displacememt Calculation  |'
write(*,*) '------------------------------------------'
write(*,*) 'ION file : ', trim(filename)//'.ion'
write(*,*) 'MSD file : ', trim(filename)//'.msd'
!
!
!Initiation
!amass(1)=15.9994_DP      ! mass of Oxygen (amu) 
!amass(2)=1.00794_DP      ! mass of H (amu)
!************************************
  
end subroutine init
!******************************************************************************************
!
!******************************************************************************************
subroutine files1
        !
        implicit none
        !
        open(unit=1, file=(trim(filename)//'.ion'), status='old')
        !open(unit=2, file=(trim(filename)//'.error'), status='old')
        open(unit=3, file=(trim(filename)//'.msd'), status='unknown')
        open(unit=4, file=(trim(filename)//'.msd_error'), status='unknown')
        !open(unit=5, file=(trim(filename)//'.test'), status='unknown')
        !
end subroutine files1

subroutine files2
        implicit none
        close(1)
        close(2)
        close(3)
        close(4)
        !close(5)
end subroutine files2 
!
!******************************************************************************************
subroutine main
!
implicit none
!
!************************************
!
!
!read from .fss to get the total ncount
ncount=0
do 
  !*********just want to know how many steps there is******************
  read(1,*,iostat=ierror) step, dummyR
  if (ierror .lt. 0) then
  !  write(*,*) '  End of File Reached'
  !  write(*,fmt='(1X, "  Total number of Samples: ", I7)' )(ncount)
    exit
  endif
  ncount=ncount+1
enddo
close(1)
open(unit=1, file=(trim(filename)//'.ion'), status='old')
!
allocate(time(ncount-stepstart+1),rIO(3,ncount-stepstart+1),box(3,ncount-stepstart+1))
allocate(drp(3,ncount-stepstart+1),step(ncount-stepstart+1))
!***********************
!
!
!
rIO=0
time=0
box=0
drp=0
step=0
e_step=0
!
!
!
!***********************
!Read .error file
!ncount=0
!read(2,*)
!do 
!  ncount=ncount+1
!  read(2,*,iostat=ierror) e_step(ncount), dummyR, e_iiO(ncount)
!  if (ierror .lt. 0) exit
!enddo
!
!
! 
!***********************
!Read .ion file
!
!********
!Jump out the first few lines
if (stepstart .gt. 1) then
  do k=1, stepstart-1
    do i=1,3
      read(1,*)
    enddo
  enddo
endif
!********
!
ncount=0
k0=1
do 
  if (cs .eq. 2) &
    read(1,*,iostat=ierror) step(ncount+1), time(ncount+1), numION, iiO, dummyI, rIO(1:3,ncount)
  if (cs .eq. 3) &
    read(1,*,iostat=ierror) step(ncount+1), time(ncount+1), numION, iiO, dummyI, dummyI, dummyI,&
                            rIO(1:3,ncount)
  !
  if (ierror .lt. 0) then
    write(*,*)
    write(*,*) '  End of File Reached'
    write(*,fmt='(1X, "  Total number of Samples: ", I7)' )(ncount)
    exit
  endif
  !
  ncount=ncount+1
  !
  do j=1,3
    if (rIO(j,ncount) .gt. 0) rIO(j,ncount)=rIO(j,ncount)-int(rIO(j,ncount)/celldm)*celldm
    if (rIO(j,ncount) .lt. 0) rIO(j,ncount)=rIO(j,ncount)-(int(rIO(j,ncount)/celldm)-1)*celldm
  enddo
  !
  if ((step(ncount)-step(ncount-1)) .ne. dstep) then
    write(*,'("Find dstep change at step=",I7,"from",I6,"to",I6, "at time", F8.4, "step", I7)') & 
         ncount, dstep,(step(ncount)-step(ncount-1)), time(ncount), step(ncount)
    dstep=step(ncount)-step(ncount-1)
    k0=ncount
  endif
  !
enddo
!
write(*,'(3x,"ncount=",I7,5x,"k0=",I3)') ncount,k0
write(*,*)
!
!
!
!
write(*,'("MSD is calculated during time of:",2F8.3)') time(k0),time(ncount)
!
!
!
!
!
!************CORE CONTENTS***************
!Start MSD calculation
!do dk=1,ncount-k0
do dk=1,last_dk
  !
  n=0
  average_box=0
  tot_sd=0
  !
  do k=k0,ncount-dk,5
    !
    do j=1,3
      !
      tmp_dr=rIO(j,k+dk)-rIO(j,k)
      dr(j)=tmp_dr+box(j,k)*celldm
      !dr=dr-nint(dr/celldm)*celldm
      !
      !
      !Adjust the parameter: box
      if (abs(dr(j)-drp(j,k)) .gt. cutoff1 ) then
        if (dr(j)-drp(j,k) .gt. 0) then
          box(j,k)=box(j,k)-1
        else 
          box(j,k)=box(j,k)+1
        endif
        dr(j)=tmp_dr+box(j,k)*celldm
      endif
      !
      !
      !Situations that can go wrong
      if (abs(dr(j)-drp(j,k)) .gt. cutoff2 ) then
        !if (                                &
        !    step(k+dk) .ne. e_step(1) .and. &
        !    step(k+dk) .ne. e_step(2) .and. &
        !    step(k+dk) .ne. e_step(3) .and. &
        !    step(k+dk) .ne. e_step(4) .and. &
        !    step(k+dk) .ne. e_step(5)       &
        !   ) then
        !  !Error! dr-drp too big!
          write(4,*) step(k+dk),time(k+dk),abs(dr(j)-drp(j,k)),j
        !endif
      endif
      !
      !
      average_box(j)=average_box(j)+box(j,k)
      drp(j,k)=dr(j)
      !
    enddo
    !
    d2=(dr(1))**2+(dr(2))**2+(dr(3))**2
    tot_sd=tot_sd+d2
    n=n+1
    !
    !
    !
    !##### Checking Errors #####
    !if (mod(dk,1000) .eq. 0) write(5,'(2I8,5x,3F10.5,5x,3I3)') k, dk, rIO(1:3,k+dk),box(1:3,k)
    !##### Checking Errors Completes#####
    !
    !
    !
  enddo
  !
  msd=tot_sd*convertBA*convertBA/real(n)
  write(3,'(2F8.3)') time(dk)-time(k0), msd
  !
enddo
!
!
!
write(*,*) "Program ends!"
write(*,'("Final box adjustments are:",3F8.3)') real(average_box(1))/n,real(average_box(2))/n,real(average_box(3))/n
!
!
!
end subroutine main
!
!******************************************************
end program mean_square_displacement
