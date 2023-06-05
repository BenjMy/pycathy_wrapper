C
C**************************  RELXOM ************************************
C
C  calculate relaxation parameter for the case NLRELX=2, using 
C  Huyakorn et al's adaptation (WRR 1986 22(13), pg 1795) of Cooley's 
C  empirical scheme (WRR 1983 19(5), pg 1274)
C
C***********************************************************************
C
      SUBROUTINE RELXOM(N,ITER,OMEGA,OMEGAP,PIKMXV,PNEW,POLD)
C
      IMPLICIT  NONE
      INTEGER   K
      INTEGER   N,ITER
      REAL*8    ZETA,DIFF,DIFMX,DIFAB,DIFABM,DIFMXP
      REAL*8    OMEGA,OMEGAP
      REAL*8    PIKMXV(*),PNEW(*),POLD(*)
C
      IF (ITER .EQ. 1) THEN
         OMEGA=1.0D0
      ELSE
         DIFABM=0.0D0
         DO K=1,N
            DIFF=PNEW(K)-POLD(K)
            DIFAB=DABS(DIFF)
            IF (DIFAB .GE. DIFABM) THEN
               DIFABM=DIFAB
               DIFMX=DIFF
            END IF
         END DO
         DIFMXP=PIKMXV(ITER-1)
         ZETA=DIFMX/(OMEGAP*DIFMXP)
         IF (ZETA .GE. -1.0D0) THEN
            OMEGA=(3.0D0 + ZETA)/(3.0D0 + DABS(ZETA))
         ELSE
            OMEGA=0.5D0/DABS(ZETA)
         END IF
      END IF
      OMEGAP=OMEGA
C
      RETURN
      END
