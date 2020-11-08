subroutine rdinp(nlines)
use mpi
implicit none

integer              :: ierr,kproc, nproc
double precision     :: init,endt,dt
character (LEN = 50) :: nlines
integer, parameter   :: iunit = 7

open(iunit, file='go', status="old")

  read(iunit,* )
  read(iunit,* ) nlines

close(iunit)



return
end subroutine
