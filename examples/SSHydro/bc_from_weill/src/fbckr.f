C
C**************************  FBCKR *************************************
C
C  calculate relative hydraulic conductivity using Brooks-Corey 
C  characteristic equation
C
C***********************************************************************
C
      DOUBLE PRECISION FUNCTION FBCKR(PSI)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      REAL*8   PSI
      INCLUDE 'SOILCHAR.H'
C
      IF (PSI .LT. BCPSAT) THEN
         FBCKR=abs(BCPSAT/PSI)**BC23B
      ELSE
         FBCKR=1.0D0
      END IF
C
      RETURN
      END
