C
C**************************  FHUKR2 ************************************
C
C  calculate relative hydraulic conductivity using Huyakorn 
C  characteristic equation, with Kr=Se**n relationship
C
C***********************************************************************
C
      DOUBLE PRECISION FUNCTION FHUKR2(PSI,SE)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      REAL*8   PSI,SE
      INCLUDE 'SOILCHAR.H'
C
      IF (PSI .LT. HUPSIA) THEN
         FHUKR2=SE**HUN
      ELSE
         FHUKR2=1.0D0
      END IF
C
      RETURN
      END
