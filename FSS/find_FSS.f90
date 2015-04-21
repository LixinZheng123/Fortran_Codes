program find_FSS
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
                           stepstop, readstep, e_step
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
integer                 :: ierror,ncount=0
integer                 :: h_num(4,15)           !For each Oxygen in FSS
                                                !there should be 1 or 2 hydrogen linked to it
real                    :: d,           &       !distance
                           tau_O(3,15), tau_H(3,30)   !position of each FSS atoms
real                    :: time
real,allocatable        :: rO(:,:),rH(:,:)
real                    :: r_iiO(3,5)
real                    :: per
!
namelist /input/  cs,filename, nsp1, nsp2, stepstop
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
write(*,*) 'Position file      : ', trim(filename)//'.pos'
write(*,*) 'Defect file        : ', trim(filename)//'.ion_angs'
write(*,*) '1st solvation file2: ', trim(filename)//'.fss_angs'
write(*,*) 'xsf file           : ', trim(filename)//'.fss_xsf'
write(*,*) 'Ion transfer file  : ', trim(filename)//'.trans'
write(*,*) 'Ion error file     : ', trim(filename)//'.error'
!
!allocation
allocate(rO(3,nsp1))
allocate(rH(3,nsp2))
!Initiation
per=0
iiO=0
iiO_record=0
iiO_p=0
iiO_p_record=0
iiHp=0
if (cs .eq. 2) numH=1
if (cs .eq. 3) numH=3
!
end subroutine init
!
!=========================================================================================
!=========================================================================================
!
subroutine files1
implicit none
open(unit=1, file=(trim(filename)//'.pos'), status='old')
open(unit=2, file=(trim(filename)//'.ion_angs'), status='unknown')
open(unit=11, file=(trim(filename)//'.fss_angs'), status='unknown')
open(unit=12, file=(trim(filename)//'.fss_xsf'), status='unknown')
open(unit=4, file=(trim(filename)//'.trans'), status='unknown')
open(unit=5, file=(trim(filename)//'.error'), status='unknown')
open(unit=10, file=(trim(filename)//'.trans_no_r'), status='unknown')
write(5,*) "This file is to record these occasions when ion has gone out of the first solvation shell of previous step."
end subroutine files1
!
subroutine files2
implicit none
close(1)
close(2)
close(3)
close(4)
close(5)
close(10)
close(11)
close(12)
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
!main loop
do 
  !**************************************************************
  !end of file
  !**************************************************************
  read(1,*,iostat=ierror) readstep, time 
  if (readstep .lt. readstep_p .and. readstep_p .ne. 0) then
    write(*,*) "Error! Readstep decreasing!"
    exit
  endif
  if (ierror .lt. 0) then
    write(*,*) '  End of File Reached'
    write(*,fmt='(1X, "  Total number of Samples: ", I7)' )(ncount)
    exit
  endif
  ncount=ncount+1
  !**************************************************************
  !**************************************************************
  !
  !MODULE 0: READ FILE .POS
  !
  !**************************************************************
  !**************************************************************
  do i = 1, nsp1
    !Read position of Oxygen
    read(1,*) rO(1:3,i)
  enddo
  do i = 1,nsp2
    !read position of Hydrogen
    read(1,*) rH(1:3,i)
  enddo
  !**************************************************************
  !**************************************************************
  !
  !MODULE 1: FIND THE ION
  !ONLY numION, INDEX iiO(1), iiH(1:numION,1), AND COORDINATES r_iiO(1:3,1) NEED TO BE PASSED TO NEXT MODULE.
  !
  !**************************************************************
  !**************************************************************
  numION = 0
  numW=0
  p=0
  e_step=0
  do iO=1,nsp1
    q=0
    do iH=1, nsp2
      call get_distance(rO(1:3,iO),rH(1:3,iH),d)
      if (d*convertBA .lt. r_cov) then
        q=q+1
        iiH(q,p+1)=iH
      endif
    enddo
    !*************************
    if (q .eq. 2) then
      numW=numW+1
    else if (q .eq. numH) then
      p=p+1
      iiO(p)=iO
    end if
    !**************************
  enddo
  numION=p
  do i=1, numION
    iiO_record(i)=iiO(i)
  enddo
  !**************************
  !1. The most straight forward situation
  !**************************
  if ((numION+numW) .eq. nsp1 .and. numION .eq. 1) then
    do j=1,3
      r_iiO(j,1) = rO(j,iiO(1))
    enddo
  endif
  !**************************
  !2. The different situations that would probably go wrong.
  !Just skip the judgement part, and adopt the previous ion index.
  !**************************
  if ((numION+numW) .ne. nsp1 .or. numION .gt. 2 .or. numION .lt. 1) then 
    !if ((numION+numW) .ne. nsp1 .or. numION .gt. 2) then
    !  write(*,*) time, numION, numION+numW
    !endif
    iiO(1)=iiO_p(1)
    do i=1,numH
      iiH(i,1)=iiHp(i)
    enddo
    do j=1,3
      r_iiO(j,1)=rO(j,iiO(1))
    enddo
  endif
  !**************************
  !3. When numION=2, try to find the real ion
  !**************************
  if ((numION+numW) .eq. nsp1 .and. numION .eq. 2) then
    do p=1,2
      do j=1,3
        r_iiO(j,p) = rO(j,iiO(p))
      enddo
    enddo
    !**************************
    !compare and change order
    !so that the total distances from O to H of iiO(1) is the shortest
    !If cs=2, OH-
    if (cs .eq. 2) then
      do p=1,numION
        call get_distance(r_iiO(1:3,p),rH(1:3,iiH(1,p)),dr(p))
      enddo
      do p=2,numION
        if(dr(p) .lt. dr(1)) then
         iiO(1)=iiO(p)
          dr(1)=dr(p)
          do j=1,3
            r_iiO(j,1)=r_iiO(j,p)
            iiH(1,1)=iiH(1,p)
          enddo
        endif
      enddo
    !IF cs=3, Find the mutual iiH that ion1 and ion2 has
    else if (cs .eq. 3) then
      do q=1,3
        do q1=1,3
          if (iiH(q,1) .eq. iiH(q1,2)) then
            iiH_real=iiH(q,1)
            exit
          endif
        enddo
      enddo
      do p=1,2
        call get_distance(r_iiO(1:3,p),rH(1:3,iiH_real),dr(p))
      enddo
      do p=2,numION
        if(dr(p) .lt. dr(1)) then
         iiO(1)=iiO(p)
          dr(1)=dr(p)
          do j=1,3
            r_iiO(j,1)=r_iiO(j,p)
            do q=1,numH
              iiH(q,1)=iiH(q,p)
            enddo
          enddo
        endif
      enddo
    endif
  endif
  !
  !
  !***
  !Modified 20140901
  !To determine if the ion has transfered too fast
  error=0
  !
  if (ncount .gt. 1) then
    do i=1,new_nsp1
      if (iiO(1) .eq. iFSS_O(i) ) then
        error=1
        exit
      endif
    enddo
  endif
  if (error .eq. 0 .and. ncount .gt. 1) then
    !write(*,*) "Error! The new ion has gone out of the FSS of previous ion at readstep=",readstep
    e_step=readstep
    write(5,*) readstep,time,iiO(1) 
  endif
  !***
  !
  !
  !**************************
  !Write in .ion
  !**************************
  !Modified 20150107 format
  ! The 2 lines below have been moved to end of loop.
  !if (cs .eq. 2) write(2,fmt='(I8,F9.4,3I5,3F9.4)' ) &
  !               readstep,time,numION,iiO(1),iiH(1:numH,1),r_iiO(1:3,1)
  if (cs .eq. 3) write(2,fmt='(I8,F9.4,5I5,3F9.4)' ) &
                 readstep,time,numION,iiO(1),iiH(1:numH,1),r_iiO(1:3,1)*convertBA
  !
  !**************************
  !**************************
  if (iiO(1) .ne. iiO_p(1)) then
    !Write in .trans
    !Modified 20150107 format
    if (cs .eq. 2) write(4,fmt='(I8,F9.4,2I5)' ) &
                   readstep, time, iiO(1), iiH(1:numH,1)
    if (cs .eq. 3) write(4,fmt='(I8,F9.4,4I5)' ) &
                   readstep, time, iiO(1), iiH(1:numH,1)
    !Modified 20141014
    !Modified 20150107 (Wrong variable name previously)
    !We also want to write out the non-rattling-transfer index in the same time.
    if (iiO_p(2) .ne. iiO(1)) then
      if (iiO_p_record(1) .ne. iiO(1) .and. iiO_p_record(2) .ne. iiO(1)) then
        write(10,fmt='(I8,I5)') readstep, iiO(1)
        iiO_p_record(2)=iiO_p_record(1)
        iiO_p_record(1)=iiO(1)
      endif
    endif
    !End of modification
    iiO_p(2)=iiO(1) 
    iiO_p(1)=iiO(1)
  endif
  do i=1,numH
    iiHp(i)=iiH(i,1)
  enddo
  !
  !
  !
  !**************************************************************
  !**************************************************************
  !
  !MODULE 2: FIND THE 1ST SOLVATION SHELL 
  !ONLY numION, INDEX iiO(1), iiH(1:numION,1), AND COORDINATES r_iiO(1:3,1) NEED TO BE PASSED FROM LAST MODULE.
  !
  !**************************************************************
  !**************************************************************
  !
  !**************************
  !This if statement is to make sure that the 1st ion has been found;
  !If not, we will write nothing to the .fss file.
  !**************************
  if (start .eq. 0 .and. numION .eq. 0) then
    cycle
  else if (start .eq. 0 .and. numION .ne. 0) then
    start=1
  endif
  !**************************
  !
  !**************************
  !Find FSS oxygens
  !**************************
  p=0
  q=0
  h_num=0
  do i=1,nsp1
    call get_distance(rO(1:3,i),r_iiO(1:3,1),d)
    if (d*convertBA .lt. r_FSS) then
      p=p+1  
      iFSS_O(p)=i
      do j=1,3
        tau_O(j,p)=rO(j,i)
      enddo
    endif
  enddo
  new_nsp1=p
  !**************************
  !Find the hydrogens relates to FSS oxygens
  !**************************
  do i=1,nsp2
    call get_distance(rH(1:3,i),r_iiO(1:3,1),dr(1))
    do p=1,new_nsp1
      call get_distance(rH(1:3,i),tau_O(1:3,p),dr(2))
      if (dr(2)*convertBA .lt. r_cov .or. dr(1)*convertBA .lt. (r_FSS-r_cov)) then
        q=q+1
        iFSS_H(q)=i
        do j=1,3
          tau_H(j,q)=rH(j,i)
        enddo
        exit
      endif     
    enddo
  enddo
  new_nsp2=q
  !**************************
  !Record the link between O and H
  !**************************
  do p=1,new_nsp1
    k=0
    do q=1, new_nsp2 
      call get_distance(tau_H(1:3,q),tau_O(1:3,p),d)
      if (d*convertBA .lt. r_cov) then
        k=k+1
        h_num(k,p)=iFSS_H(q)
      endif
    enddo
  enddo
  !**************************
  !Write index into file.fss 
  !**************************
  if (numION .eq. 1) write(11,*) readstep, time, numION, new_nsp1,new_nsp2
  if (numION .ne. 1) write(11,*) readstep, time, numION, new_nsp1,new_nsp2, iiO_record(1:numION)
  write(12,*) readstep, time
  write(12,*) new_nsp1+new_nsp2, "      1"
  write(11,*) iiO(1),iiH(1:numH,1)
  do p=1,new_nsp1
    do j=1,3
      if (tau_O(j,p) .gt. celldm) tau_O(j,p)=tau_O(j,p)-int(tau_O(j,p)/celldm)*celldm
      if (tau_O(j,p) .lt. 0) tau_O(j,p)=tau_O(j,p)+(abs(int(tau_O(j,p)/celldm))+1)*celldm
    enddo
    write(11,*) iFSS_O(p),tau_O(1:3,p)*convertBA,h_num(1:4,p)
    write(12,*) "O   ",tau_O(1:3,p)*convertBA
  enddo
  do q=1,new_nsp2
    do j=1,3
      if (tau_H(j,q) .gt. celldm) tau_H(j,q)=tau_H(j,q)-int(tau_H(j,q)/celldm)*celldm
      if (tau_H(j,q) .lt. 0) tau_H(j,q)=tau_H(j,q)+(abs(int(tau_H(j,q)/celldm))+1)*celldm
    enddo
    write(11,*) iFSS_H(q),tau_H(1:3,q)*convertBA
    write(12,*) "H   ",tau_H(1:3,q)*convertBA
    !Modified 20150204
    if (cs .eq. 2) then
      if(iFSS_H(q) .eq. iiH(1,1)) write(2,fmt='(I8,F9.4,3I5,6F9.4)') &
                    readstep,time,numION,iiO(1),iiH(1,1),r_iiO(1:3,1)*convertBA,tau_H(1:3,q)*convertBA
    endif
  enddo
  !
  !
  !Modified 20140901
  if (e_step .eq. readstep) then
    if (cs .eq. 2) then
      if ((new_nsp1*2-1) .ne. new_nsp2) then
        write(5,*) "Strange nsp in FSS Calculation!"
        write(5,*) readstep, time, numION, new_nsp1,new_nsp2, iiO_record(1:numION)
      else if ((new_nsp1*2+1) .ne. new_nsp2) then
        write(5,*) "Strange nsp in FSS Calculation!"
        write(5,*) readstep, time, numION, new_nsp1,new_nsp2, iiO_record(1:numION)
      endif
    endif
  endif

  readstep_p=readstep

enddo

!
end subroutine main
!
!=========================================================================================
!=========================================================================================
!
subroutine get_distance(pos1,pos2,dist)

implicit none
!
real, intent(in)    :: pos1(3),pos2(3)
real, intent(inout) :: dist
real                :: delta(3)
real, parameter     :: pcell=23.5170
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
!
!=========================================================================================
!=========================================================================================
!
end program find_FSS
