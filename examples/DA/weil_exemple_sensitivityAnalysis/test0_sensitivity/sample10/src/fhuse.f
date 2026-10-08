C
C**************************  FHUSE *************************************
C
C  calculate effective saturation using Huyakorn characteristic equation
C
C***********************************************************************
C
      DOUBLE PRECISION FUNCTION FHUSE(PSI)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      REAL*8   LAMBDA
      REAL*8   PSI
      INCLUDE 'SOILCHAR.H'
C
      IF (PSI .LT. HUPSIA) THEN
         LAMBDA=HUALB*((HUPSIA - PSI)**HUBETA)
         FHUSE=(1.0D0/(1.0D0 + LAMBDA))**HUGAMA
      ELSE
         FHUSE=1.0D0
      END IF
C
      RETURN
      END
