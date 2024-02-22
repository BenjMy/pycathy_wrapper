C
C**************************  FHUDSE ************************************
C
C  calculate derivative of effective saturation with respect to pressure 
C  head using Huyakorn characteristic equation
C
C***********************************************************************
C
      DOUBLE PRECISION FUNCTION FHUDSE(PSI)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      REAL*8   PAP,LAMBDA,LAMR
      REAL*8   PSI
      INCLUDE 'SOILCHAR.H'
C
      IF (PSI .LT. HUPSIA) THEN
         PAP=HUPSIA - PSI
         LAMBDA=HUALB*(PAP**HUBETA)
         LAMR=1.0D0/(1.0D0 + LAMBDA)
         FHUDSE=(HUGB*LAMBDA/PAP)*(LAMR**HUGAM1)
      ELSE
         FHUDSE=0.0D0
      END IF
C
      RETURN
      END
