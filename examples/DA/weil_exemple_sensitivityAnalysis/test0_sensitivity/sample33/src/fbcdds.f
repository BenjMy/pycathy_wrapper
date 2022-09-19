C
C**************************  FBCDDS ************************************
C
C  calculate second derivative of effective saturation with respect to
C  pressure head using Brooks-Corey characteristic equation
C
C***********************************************************************
C
      DOUBLE PRECISION FUNCTION FBCDDS(PSI)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      REAL*8   PSI
      INCLUDE 'SOILCHAR.H'
C
      IF (PSI .LT. BCPSAT) THEN
         FBCDDS=BCBB1P*(BCPSAT/PSI)**BCB2
      ELSE
         FBCDDS=0.0D0
      END IF
C
      RETURN
      END
