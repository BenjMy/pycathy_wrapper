C
C**************************  RELAX  ************************************
C
C  apply nonlinear relaxation parameter OMEGA to update current 
C  pressure head values
C
C***********************************************************************
C
      SUBROUTINE RELAX(N,OMEGA,PNEW,POLD)
C
      IMPLICIT  NONE
      INTEGER   K
      INTEGER   N
      REAL*8    OMEGA1
      REAL*8    OMEGA
      REAL*8    PNEW(*),POLD(*)
C
      OMEGA1=1.0D0 - OMEGA
      DO K=1,N
         PNEW(K)=OMEGA*PNEW(K) + OMEGA1*POLD(K)
      END DO
C
      RETURN
      END
