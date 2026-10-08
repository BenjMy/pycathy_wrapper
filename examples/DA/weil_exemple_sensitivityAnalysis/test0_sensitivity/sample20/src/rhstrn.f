C
C**************************  RHSTRN ************************************
C
C  assemble RHS vector, without the contributions of the 
C  boundary conditions
C
C***********************************************************************
C
      SUBROUTINE RHSTRN(N,TOPOLC,JAC,DELTAT,TETAC,TNOTIC,CTIMEP,
     1                  COEF1C,COEF2C)
C
      IMPLICIT  NONE
      INTEGER   I,J,K,M
      INTEGER   N
      INTEGER   TOPOLC(*),JAC(*)
      REAL*8    RDT,TETA1,COEF
      REAL*8    DELTAT,TETAC
      REAL*8    TNOTIC(*),CTIMEP(*),COEF1C(*),COEF2C(*)
C
      call init0r(n,tnotic)
      RDT=1.0D0/DELTAT
      TETA1=1.0D0 - TETAC
      DO K=1,N
        I=TOPOLC(K)
        J=TOPOLC(K+1)-1
        DO M=I,J
          COEF=COEF2C(M)*RDT - TETA1*COEF1C(M)
          TNOTIC(K)=TNOTIC(K) + COEF*CTIMEP(JAC(M))
          IF(M.NE.I) TNOTIC(JAC(M))=TNOTIC(JAC(M))+COEF*CTIMEP(K)
        END DO
      END DO


C
      RETURN
      END
