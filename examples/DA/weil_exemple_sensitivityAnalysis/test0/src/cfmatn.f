C
C**************************  CFMATN ************************************
C
C  assemble the global LHS system matrix (the Jacobian) from the 
C  stiffness and mass matrices and the derivative term components
C  of the Jacobian: Newton scheme
C
C***********************************************************************
C
      SUBROUTINE CFMATN(NTERM,TETAF,DELTAT,COEF1,COEF2,COEF3)
C
      IMPLICIT  NONE
      INTEGER   K
      INTEGER   NTERM
      REAL*8    RDT
      REAL*8    TETAF,DELTAT
      REAL*8    COEF1(*),COEF2(*),COEF3(*)
C
      RDT=1.0D0/DELTAT
      DO K=1,NTERM
         COEF1(K)=TETAF*COEF1(K) + COEF2(K)*RDT + COEF3(K)
      END DO
C
      RETURN
      END
