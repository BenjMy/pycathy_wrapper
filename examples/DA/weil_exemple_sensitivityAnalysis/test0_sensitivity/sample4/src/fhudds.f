C
C**************************  FHUDDS ************************************
C
C  calculate second derivative of effective saturation with respect to
C  pressure head using Huyakorn characteristic equation
C
C***********************************************************************
C
      DOUBLE PRECISION FUNCTION FHUDDS(PSI)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      REAL*8   PAP,PAPR,LAMBDA,LAMR
      REAL*8   PSI
      INCLUDE 'SOILCHAR.H'
C
      IF (PSI .LT. HUPSIA) THEN
         PAP=HUPSIA - PSI
         PAPR=1.0D0/PAP
         LAMBDA=HUALB*(PAP**HUBETA)
         LAMR=1.0D0/(1.0D0 + LAMBDA)
         FHUDDS=HUGB*LAMBDA*PAPR*PAPR*(HU1BET + HUGB1*LAMBDA)*
     1                                (LAMR**HUGAM2)
      ELSE
         FHUDDS=0.0D0
      END IF
C
      RETURN
      END
