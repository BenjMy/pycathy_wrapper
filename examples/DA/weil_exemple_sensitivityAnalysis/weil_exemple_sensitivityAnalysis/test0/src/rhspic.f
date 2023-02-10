C
C**************************  RHSPIC ************************************
C
C  assemble RHS vector, without contribution of the unsaturated
C  zone gravitational term and without the boundary conditions:
C  Picard scheme
C
C***********************************************************************
C
      SUBROUTINE RHSPIC(N,JA,TOPOL,DELTAT,PTNEW,PNEW,PTIMEP,
     1                 SWNEW,SWTIMEP,COEF1,COEF2,COEF4,TNOTI)
C
      IMPLICIT  NONE
      INTEGER   I,J,K,M
      INTEGER   N
      INTEGER   JA(*),TOPOL(*)
      REAL*8    RDT,COEF,COEFF
      REAL*8    DELTAT
      REAL*8    PTNEW(*),PNEW(*),PTIMEP(*),COEF1(*),COEF2(*),TNOTI(*)
      REAL*8    SWNEW(*),SWTIMEP(*),COEF4(*)
C
      CALL INIT0R(N,TNOTI)
      RDT=1.0D0/DELTAT
      DO K=1,N
         I=TOPOL(K)
         J=TOPOL(K+1)-1
         DO M=I,J
            COEF=COEF2(M)*RDT 
            COEFF=COEF4(M)*RDT
            TNOTI(K)=TNOTI(K) - COEF1(M)*PTNEW(JA(M)) - 
     1                          COEF*(PNEW(JA(M)) - PTIMEP(JA(M))) -
     2                          COEFF*(SWNEW(JA(M))-SWTIMEP(JA(M)))       
            IF (M .NE. I) TNOTI(JA(M))=TNOTI(JA(M)) - 
     1                          COEF1(M)*PTNEW(K) - 
     2                          COEF*(PNEW(K) - PTIMEP(K)) -
     3                          COEFF*(SWNEW(K)-SWTIMEP(K))       
         END DO
      END DO
C
      RETURN
      END
