C
C**************************  FHUKR3 ************************************
C
C  calculate relative hydraulic conductivity using Huyakorn 
C  characteristic equation, with conductivity relationship from 
C  Table 3 of 1984 paper (log_10 Kr(Se) curve)
C
C***********************************************************************
C
      DOUBLE PRECISION FUNCTION FHUKR3(PSI,SE)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      REAL*8   PSI,SE
      INCLUDE 'SOILCHAR.H'
C
      IF (PSI .LT. HUPSIA) THEN
         FHUKR3=10.0D0**(HUA*SE*SE + HUB2A*SE + HUAB)
      ELSE
         FHUKR3=1.0D0
      END IF
C
      RETURN
      END
