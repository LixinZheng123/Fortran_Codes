!Radial Distribution Function g(r) for Quantum Espresso Output Filesres
!Calculates wannier function center distribution around a certain O*.
!Lixin Zheng, October 2014
!==================================================================================================
program rdf 
   !
   implicit none
   !
   INTEGER, PARAMETER   :: DP = selected_real_kind(14,200) !double-precision kind
   real(DP), parameter  :: pi=4.d0*atan(1.d0),&
                           celldm=23.5170,    &
                           convertBA = 0.52917
   character(len=1)     :: rs
   integer              :: binNum,        &  !number of bins for g(r)
                           nsp(2),        &  !number of each atomic species     
                           cs,            &  !cs=2, OH-; cs=3, H3O+
                           stepstart=1,   &  !nfi that we will start reading *.pos 
                           stepstop=-1,   &  !nfi that we will stop reading *.pos AFTER (optional)
                           ncount=0,      &  !number of steps index
                           i, j, k, ii=0, &  !general indexes
                           n,m,l,p,       &  !general index 
                           step, stepp,   &  !read the nfi from *.pos 
                           nbin              !bin index
   integer              :: ierror
   real(DP), allocatable:: grt(:),        &  !unnormalized g(r)
                           temgrt(:)         !temp unnormalized g(r)
   real(DP)             :: gofr,          &  !FINAL g(r) print-out variable (not an array, print-out values)  
                           sum=0,         &
                           r,             &  !total distance (not an array, print-out value)
                           d,             &
                           cell(3,3),     &  !lattice prim vectors
                           cellinv(3,3),  &  !inverse of prim vectors, for scaled positions
                           rO(3,64),rw(3,256),&!atomic positions (dimension,atomic-species)
                           box,           &  !length of the region (includes multiple supercells)
                           res,           &  !resolution (size) of bins  binNum/Hbox (number/length)
                           omega,         &  !volume of cell
                           dens,          &  !raw particle density
                           r2,            &  !total squared distance
                           time,          &  !timestep
                           sigma,         &
                           vshell            !volume of the infinitesimal radial shell

   character(len=30)       :: filename, filegofr
   !
   namelist /input/ filename,filegofr, nsp, cs, binNum, stepstart, stepstop
   !
   !initialization
   call init
   !
   !Open files
   call files(1)
   !
   !Main Loop
   call main
   !
   !finalize the results
   call final
   !
   !close files
   call files(2)
   !
   !
   !********************************************************************************
   contains
   !********************************************************************************
      !
      !
subroutine init
!
implicit none
!
!read Namelist input for stdin
read(*,input)
!
!
!
!Print Intro
write(*,*) '------------------------------'
write(*,*) '|          g(r) of wfc        |'
write(*,*) '------------------------------'
write(*,*) 'Input file : ', trim(filename)//'.pos', trim(filename)//'.wfc'
write(*,*) 'Output file        : ', trim(filegofr)
write(*,fmt='(1X, "Bin number          : ", I8)') binNum
!
!Allocate  Files
allocate(grt(binNum), temgrt(binNum))
!
grt=0
!******************************************************
!calculate the inverse
!cell=0
!cell(1,1)=celldm
!cell(2,2)=celldm
!cell(3,3)=celldm
!call invert(cell, cellinv, omega)
!box=omega**(1.0d0/3.0d0)
box=1.2
omega=box*box*box
res=binNum/box
!******************************************************
end subroutine init
!
!******************************************************
subroutine files(io)
!
implicit none
!
integer, intent(in)  ::io
!
select case(io)
  case(1)
    open(unit=1, file=(trim(filename)//'.pos'), status='old')
    open(unit=2, file=(trim(filename)//'.wfc'), status='old')
    open(unit=3, file=(trim(filegofr)), status='unknown')
    !
  case(2)
    close(1)
    close(2)
    close(3)
  case default
end select
!
end subroutine files
!
!
!
!******************************************************
!------------------------------------------------------
!Main Loop
! read-in data, calculate inverse, r_to_s,
! 
! 
!------------------------------------------------------
!******************************************************
subroutine main
!
implicit none
!
!Jump out the first few lines
!if (stepstart .gt. 1) then
!  do k=1, stepstart-1
!    read(1,*) step, time
!    do i=1,nsp(1)+nsp(2)+1
!      read(1,*)
!    enddo
!  enddo
!endif
!
do 
  !******************************************************
  read(1,*,iostat=ierror) step, time
  !if (ncount .gt. stepstop .and. stepstop .ne. -1) exit
  if (ierror .lt. 0) then
    write(*,*) '  End of File Reached'
    write(*,fmt='(1X, "  Total number of Samples: ", I7)' )(ncount)
    exit
  endif
  !
  read(2,*,iostat=ierror) stepp, time
  if (ierror .lt. 0) then
    write(*,*) '  End of File Reached'
    write(*,fmt='(1X, "  Total number of Samples: ", I7)' )(ncount)
    exit
  endif
  !
  !if (step .ne. stepp) write(*,*) "ERROR! Readstep not equal in 2 files!", step,stepp
  !******************************************************
  !Read in .pos
  do i=1, nsp(1)
    read(1,*)  rO(1:3,i)
  enddo
  do i=1,nsp(2)
    read(1,*)
  enddo
  !Read in .wfc
  do i=1, nsp(1)*4
    read(2,*)  rw(1:3,i)
  enddo
  !******************************************************
  !Zero all temp g(r)
  temgrt = 0
  !******************************************************
  !------------------Count PAIRS-----------------------
  do i=1,nsp(1)
    do k=1,nsp(1)*4
      call get_distance(rO(1:3,i),rw(1:3,k),d)
      if (d .lt. box) then
        nbin = nint(d*res)
        if (nbin .gt. binNum) nbin = binNum
        temgrt(nbin) = temgrt(nbin) + 1
      endif
    enddo
  enddo
  !=====================================================
  do n=1,binNum,1
    grt(n)=grt(n)+temgrt(n)
  enddo
  ncount=ncount + 1
enddo 
!
write(*,*) ' ... Main loop completed!'
!
end subroutine main
!
!
!
!******************************************************
!------------------------------------------------------
!Finalize the results, print out to the output
!
!------------------------------------------------------
!******************************************************
subroutine final
!
implicit none
!
do n=1,binNum-1,1
  !
  !construct some details of normalization and distance  
  r = DBLE(n)/res
  vshell = 4.0*pi*(r**2)/res
  !Calculate the final, normalized, value of g(r)! Please note that the
  !normalization constant (ncount*mic3*nsp(atom1)*dens*vshell)...
  !  ncount*mic3 = number of steps and number of extended shells
  !  nsp(atom1)*dens = number of pairs
  !  vshell = volume of the infinitesimal shell
  dens=DBLE(nsp(1)*4/omega)
  gofr = grt(n)/(DBLE(ncount*vshell*nsp(1)*dens))
  sum=sum+gofr
  !
  write(3,*) (r*convertBA), gofr
  !
enddo
write(*,*) "sum of all gofr is :", sum
!
write(*,*) ' Program complete!'
!
end subroutine final
! 
!
!
!******************************************************
!******************************************************
!         
subroutine invert(Mi, Mo, det)

implicit none

real(8), intent(in)     :: Mi(3,3)
real(8), intent(out)    :: Mo(3,3), det
real(8)                 :: tmp(3,3), vs

integer              :: vi,vj,vk,vl  !indexes
integer              :: vn,vr     !int counters

det=0.0
vs   = 1.0
vi   = 1
vj   = 2
vk   = 3


do
  do vn=1,3,1
    det= det + vs*Mi(1,vi)*Mi(2,vj)*Mi(3,vk)
    vl=vi
    vi=vj
    vj=vk
    vk=vl
  enddo
  vi=2
  vj=1
  vk=3
  vs=-vs
  if (vs .GT. 0.0) exit
enddo

if(ABS(det) .LT. 1.0e-20) then
  write(*,*) 'Error: Singular Matrix'
  Stop
endif

vi=1
vj=2
vk=3

do vr=1,3
  tmp(vr,1) = (Mi(2,vj)*Mi(3,vk) - Mi(2,vk)*Mi(3,vj))/det
  tmp(vr,2) = (Mi(3,vj)*Mi(1,vk) - Mi(3,vk)*Mi(1,vj))/det
  tmp(vr,3) = (Mi(1,vj)*Mi(2,vk) - Mi(1,vk)*Mi(2,vj))/det
  vl=vi
  vi=vj
  vj=vk
  vk=vl
enddo

do vl=1,3
  do vk=1,3
    Mo(vk,vl)=tmp(vk,vl)
  enddo
enddo

end subroutine invert
!
!******************************************************
!******************************************************
!
subroutine get_distance(pos1,pos2,dist)

implicit none
!
real(DP), intent(in)    :: pos1(3),pos2(3)
real(DP), intent(inout) :: dist
real(DP)                :: delta(3)
real(DP), parameter     :: pcell=23.5170
integer                 :: i
!
do i=1,3
  delta(i)=pos1(i)-pos2(i)
  delta(i)=delta(i)-nint(delta(i)/pcell)*pcell
enddo
!
dist=sqrt(delta(1)**2+delta(2)**2+delta(3)**2)
!
return
end subroutine get_distance
!******************************************************************************************
end PROGRAM rdf 
