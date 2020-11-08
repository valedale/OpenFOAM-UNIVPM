subroutine rnum(nomef,n)
!======================================
!
implicit none
!
character (LEN=20) :: str
character (LEN=50) :: nlines
integer            :: n,ne,ok,n2,n3
character (LEN=50) :: nomef
double precision   :: n1!,nomef
integer,parameter  :: cnt=9
!
open(unit=1,file=nomef,status='old')
!
n=0
do
    read(1,*,iostat=ok) n1
    if (ok /= 0)  exit
    n = n+1
enddo
!
close(1)
!

return
end
