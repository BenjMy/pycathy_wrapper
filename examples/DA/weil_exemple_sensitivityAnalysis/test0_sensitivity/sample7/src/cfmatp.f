C
C**************************  CFMATP ************************************
C
C  assemble the global LHS system matrix from the stiffness and mass
C  matrices: Picard scheme
C
C***********************************************************************
C
      SUBROUTINE CFMATP(N,NTERM,TETAF,DELTAT,COEF1,COEF2,COEF4,COEF5,
     1                 TOPOL,ET2)
C
      IMPLICIT  NONE
      INTEGER   K
      INTEGER   NTERM,N,TOPOL(*)
      REAL*8    RDT
      REAL*8    TETAF,DELTAT
      REAL*8    COEF1(*),COEF2(*),COEF4(*),COEF5(*)
      REAL*8    ET2(*)
C
      RDT=1.0D0/DELTAT
      DO K=1,N
         COEF5(TOPOL(K))=COEF4(TOPOL(K))*ET2(K)
      END DO 
      DO K=1,NTERM
         COEF1(K)=TETAF*COEF1(K) + COEF2(K)*RDT+COEF5(K)*RDT
      END DO
C
      RETURN
      END
