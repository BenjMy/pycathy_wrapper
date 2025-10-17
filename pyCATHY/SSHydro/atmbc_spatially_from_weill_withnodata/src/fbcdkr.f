C
C**************************  FBCDKR ************************************
C
C  calculate derivative of relative hydraulic conductivity with respect
C  to pressure head using Brooks-Corey characteristic equation
C
C***********************************************************************
C
      DOUBLE PRECISION FUNCTION FBCDKR(PSI)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      REAL*8   PSI
      INCLUDE 'SOILCHAR.H'
C
      IF (PSI .LT. BCPSAT) THEN
         FBCDKR=BC23BP*(BCPSAT/PSI)**BC33B
      ELSE
         FBCDKR=0.0D0
      END IF
C
      RETURN
      END
