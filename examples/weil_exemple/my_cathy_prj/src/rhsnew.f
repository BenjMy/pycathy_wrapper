C
C**************************  RHSNEW ************************************
C
C  assemble RHS vector, without contribution of the unsaturated
C  zone gravitational term and without the boundary conditions:
C  Newton scheme
C
C***********************************************************************
C
      SUBROUTINE RHSNEW(N,JA,TOPOL,DELTAT,PTNEW,PNEW,
     1                  PTIMEP,COEF1,COEF2,TNOTI)
C
      IMPLICIT  NONE
      INTEGER   I,J,K,M
      INTEGER   N
      INTEGER   JA(*),TOPOL(*)
      REAL*8    RDT
      REAL*8    DELTAT
      REAL*8    PTNEW(*),PNEW(*),PTIMEP(*),COEF1(*),COEF2(*),TNOTI(*)
C
      RDT=1.0D0/DELTAT
      DO K=1,N
         TNOTI(K)=0.0D0
         I=TOPOL(K)
         J=TOPOL(K+1)-1
         DO M=I,J
            TNOTI(K)=TNOTI(K) - COEF1(M)*PTNEW(JA(M)) - COEF2(M)*
     1               (PNEW(JA(M)) - PTIMEP(JA(M)))*RDT
         END DO
      END DO
C
      RETURN
      END
