C
C**************************  CFMAT  ************************************
C
C  assemble the global LHS system matrix from the stiffness and mass
C  matrices 
C
C***********************************************************************
C
      SUBROUTINE CFMAT(NTERMC,TETAC,DELTAT,COEF1C,COEF2C)
C
      IMPLICIT  NONE
      INTEGER   K
      INTEGER   NTERMC,N
      REAL*8    RDT
      REAL*8    TETAC,DELTAT
      REAL*8    COEF1C(*),COEF2C(*)
C
      RDT=1.0D0/DELTAT
      DO K=1,NTERMC
         COEF1C(K)=+TETAC*COEF1C(K) + COEF2C(K)*RDT
      END DO
C
      RETURN
      END
