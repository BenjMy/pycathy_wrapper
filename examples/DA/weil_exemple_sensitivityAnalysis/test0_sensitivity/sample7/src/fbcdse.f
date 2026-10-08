C
C**************************  FBCDSE ************************************
C
C  calculate derivative of effective saturation with respect to pressure 
C  head using Brooks-Corey characteristic equation
C
C***********************************************************************
C
      DOUBLE PRECISION FUNCTION FBCDSE(PSI)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      REAL*8   PSI
      INCLUDE 'SOILCHAR.H'
C
      IF (PSI .LT. BCPSAT) THEN
         FBCDSE=BCBPS*abs(BCPSAT/PSI)**BCB1
      ELSE
         FBCDSE=0.0D0
      END IF
C
      RETURN
      END
