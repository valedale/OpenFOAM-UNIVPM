program pdft
use mpi
implicit none

integer                        :: i,j,id
integer                        :: ierr, nproc, kproc
integer                        :: nr,nrp,its
character (LEN = 50)           :: nlines
integer, parameter             :: kunit = 20
double precision, allocatable  :: dat(:,:), datp(:,:)
integer, allocatable :: nntp(:),offset(:), nrp_p(:), ioffset_p(:)
integer                        :: ioffset 
!
double precision   :: times,start, ctime, ctime_g
character (LEN=20) :: str
!-----------------------------------------------------------
!
call MPI_INIT(ierr)

start= MPI_Wtime() ; 

call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,kproc,ierr)

if ( kproc == 0 ) then
    call rdinp(nlines)
    call rnum (nlines,nr)
endif
call MPI_BCAST (nr,nproc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

allocate (dat(nr,2))
allocate (offset(nproc))
allocate (nntp(nproc))
allocate (nrp_p(nproc))
allocate (ioffset_p(nproc))

nrp  = int(nr/nproc)
its = nr  - nrp*nproc !! Time -- step disavanz

nntp(1:nproc)  =0

if ( its > 0 ) then
  if (kproc +1 <= its) then
         nrp = nrp + 1
  endif
endif

!_---------------------------------------------------------
nrp_p(:) = 0 
call MPI_Gather (nrp,1 , MPI_INTEGER, nrp_p, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr ) 

if (kproc == 0) then
ioffset_p(:) = 0
    do i = 2, nproc
     ioffset_p(i) = 0 
      do j=2, i
         ioffset_p(i) = ioffset_p(i) + nrp_p(j-1)
      enddo
    enddo
endif

call MPI_Scatter (ioffset_p, 1 , MPI_INTEGER, ioffset, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr ) 

!
!-----------------------------------------------------------
!
if (kproc .eq. 0 ) then
   open(kunit, file=nlines(1:LEN_TRIM(nlines)), status="old")

   do i= 1,nr
      read (kunit,*) dat(i,:)
   enddo
!
endif
allocate(datp(nrp,2))
call MPI_BCAST (dat,nr*2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!
!-----------------------------------------------------------
!
do id= 1,nrp
     datp(id,:) = dat(id + ioffset,:) 
enddo
deallocate(dat)
!
!-----------------------------------------------------------
!
call dftc(nr,nrp,ioffset,datp)
!
!-----------------------------------------------------------
!
write(*,*) 'my rank ',kproc, ' ioffset ',ioffset

times =  MPI_Wtime() ; 
ctime  = times  - start
call MPI_Reduce (ctime,ctime_g,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!
!-----------------------------------------------------------
!
if (kproc == 0 ) then
   write (*,16) ctime_g/nproc
endif 

call MPI_FINALIZE(ierr)
STOP

16 format(/,' clock time = ',e12.6,' s')

end program
