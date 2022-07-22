C
C**************************  FXVKR *************************************
C
C  calculate relative hydraulic conductivity using extended 
C  van Genuchten characteristic equation. FXVKR is the same function
C  as FVGKR, but computed differently.
C
C***********************************************************************
C
      DOUBLE PRECISION FUNCTION FXVKR(PSI,INOD)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  INOD
      REAL*8   BETA,B1,V1
      REAL*8   PSI
      INCLUDE 'SOILCHAR.H'
C
      IF (PSI .LT. -1.0D-14) THEN
         BETA=(PSI/VGPSAT(INOD))**VGN(INOD)
         B1=BETA+1.0D0
         V1=(B1**VGM(INOD)) - (BETA**VGM(INOD))
         FXVKR=((1.0D0/B1)**VGM52(INOD))*V1*V1
      ELSE
         FXVKR=1.0D0
      END IF
C
      RETURN
      END
