program pt
call precisiontest
end program

subroutine precisiontest()

use mpar
real(kind=RSP)   :: real_sp
real(kind=REP) :: real_ep
integer(kind=ISP):: integer_sp
integer(kind=IEP):: integer_ep

WRITE(*,*) 
WRITE(*,*) "REAL(KIND=RSP), REAL SINGLE PRECISION"
WRITE(*,*) "RSP =", RSP
WRITE(*,*) "Number of decimal digits of precision =", precision(real_sp)
WRITE(*,*) "Maximum binary exponent =", maxexponent(real_sp)
WRITE(*,*) "Minimim binary exponent =", minexponent(real_sp)
WRITE(*,*)
WRITE(*,*) "REAL(KIND=REP), REAL EXTENDED PRECISION"
WRITE(*,*) "REP =", REP
WRITE(*,*) "Number of decimal digits of precision =", precision(real_ep)
WRITE(*,*) "Maximum binary exponent =", maxexponent(real_ep)
WRITE(*,*) "Minimim binary exponent =", minexponent(real_ep)
WRITE(*,*)
WRITE(*,*) "INTEGER(KIND=ISP), INTEGER SINGLE PRECISION"
WRITE(*,*) "ISP =", ISP
WRITE(*,*) "Range =", range(integer_sp)
WRITE(*,*)
WRITE(*,*) "INTEGER(KIND=IEP), INTEGER EXTENDED PRECISION"
WRITE(*,*) "IEP =", IEP
WRITE(*,*) "Range =", range(integer_ep)
WRITE(*,*)

end subroutine precisiontest
