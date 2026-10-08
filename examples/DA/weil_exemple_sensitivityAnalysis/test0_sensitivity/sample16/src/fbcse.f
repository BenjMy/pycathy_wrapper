C
C**************************  FBCSE *************************************
C
C  calculate effective saturation using Brooks-Corey 
C  characteristic equation
C
C***********************************************************************
C
      DOUBLE PRECISION FUNCTION FBCSE(PSI)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      REAL*8   PSI
      INCLUDE 'SOILCHAR.H'
C
      IF (PSI .LT. BCPSAT) THEN
         FBCSE=abs(BCPSAT/PSI)**BCBETA
      ELSE
         FBCSE=1.0D0
      END IF
C
      RETURN
      END
