subroutine dftc(nr,nrp,ioffset,datp)
use mpi
implicit none

integer :: i,nrp,nr,nf,ioffset
integer :: j,kproc,ierr
double precision   :: datp(nrp,2),y(nrp),dt
character (LEN=20) :: str
integer, parameter :: kunit = 10
double precision,allocatable :: vec_f(:),cdft(:)
double precision, parameter  :: pi = 3.1415926535d0

complex(kind=8)  :: vec_e(nrp)
complex(kind=8),allocatable :: h(:)
complex(kind=8),allocatable :: ht(:)
!
!-----------------------------------------------------------------------------------
!
call MPI_COMM_RANK(MPI_COMM_WORLD,kproc,ierr)
!
!-----------------------------------------------------------------------------------
!
nf = (nr-1)/2
allocate ( h(nf+1))
allocate ( ht(nf+1) )
allocate ( vec_f((nr/2-1/2+1)) )
allocate ( cdft((nr/2-1/2+1)) )

do i=1,nf+1
   vec_f(i) = i-1
enddo
!
do i=1,nrp
   vec_e(i) = exp((-2d0*pi/nr)*(0,1)*(i-1 + ioffset))
enddo
!
y(:) = datp(:,2)
!
do i=1,nf+1
   h(i) = (0,0)
   do j=1,nrp
      h(i) = h(i) + y(j)*vec_e(j)**(vec_f(i))
   enddo
enddo
!-----------------------------------------------------------------------------------

call MPI_AllReduce (h,ht,nr,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,ierr)

dt = datp(2,1) - datp(1,1) 

do i=1,nf+1
   vec_f(i) = vec_f(i)/(nr*dt)
   cdft(i)  = 2d0*((ht(i))*conjg(ht(i)))/(nr**2)
enddo

if (kproc == 0) then
    open(kunit+kproc, file='file_dft.dat', status="unknown")
!
    do i=1,nf+1
      write(kunit+kproc,*) vec_f(i), cdft(i)
    enddo
!
endif

return
end subroutine
