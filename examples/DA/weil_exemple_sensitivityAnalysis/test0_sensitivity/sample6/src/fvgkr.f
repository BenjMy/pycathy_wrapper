C
C**************************  FVGKR *************************************
C
C  calculate relative hydraulic conductivity using van Genuchten 
C  characteristic equation
C
C***********************************************************************
C
      DOUBLE PRECISION FUNCTION FVGKR(PSI,SE,INOD)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  INOD
      REAL*8   OMEGA,V1
      REAL*8   PSI,SE
      INCLUDE 'SOILCHAR.H'
C
      IF (PSI .LT. -1.0D-14) THEN
         OMEGA=ABS(SE)**VGMR(INOD)
         V1=1.0D0 - (ABS(1.0D0 - OMEGA)**VGM(INOD))
         FVGKR=DSQRT(SE)*V1*V1
      ELSE
         FVGKR=1.0D0
      END IF
C
      RETURN
      END
