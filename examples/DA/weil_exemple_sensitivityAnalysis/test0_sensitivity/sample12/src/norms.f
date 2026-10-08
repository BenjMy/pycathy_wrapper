C
C**************************  NORMS  ************************************
C
C  calculate the nonlinear convergence and residual errors, using
C  the L2 and infinity norms
C
C***********************************************************************
C
      SUBROUTINE NORMS(N,IKMAX,PIKMAX,PINF,PL2,FINF,FL2,PNEW,POLD,TNOTI)
C
      IMPLICIT  NONE
      INTEGER   K
      INTEGER   N,IKMAX
      REAL*8    DIFF,DIFABS,FABS
      REAL*8    PIKMAX,PINF,PL2,FINF,FL2
      REAL*8    PNEW(*),POLD(*),TNOTI(*)
C
      PL2=0.0D0
      FL2=0.0D0
      PINF=0.0D0
      FINF=0.0D0
      DO K=1,N
         DIFF=PNEW(K)-POLD(K)
         DIFABS=DABS(DIFF)
         PL2=PL2 + DIFF*DIFF
         FL2=FL2 + TNOTI(K)*TNOTI(K)
         IF (DIFABS .GE. PINF) THEN
            PINF=DIFABS
            PIKMAX=DIFF
            IKMAX=K
         END IF
         FABS=DABS(TNOTI(K))
         IF (FABS .GE. FINF) THEN
            FINF=FABS
         END IF
      END DO
      PL2=DSQRT(PL2)
      FL2=DSQRT(FL2)
C
      RETURN
      END
