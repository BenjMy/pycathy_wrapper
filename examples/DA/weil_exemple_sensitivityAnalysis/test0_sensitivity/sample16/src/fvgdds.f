C
C**************************  FVGDDS ************************************
C
C  calculate second derivative of effective saturation with respect to
C  pressure head using van Genuchten characteristic equation
C
C***********************************************************************
C
      DOUBLE PRECISION FUNCTION FVGDDS(PSI,INOD)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  INOD
      REAL*8   BETA,B1,B1R
      REAL*8   PSI
      INCLUDE 'SOILCHAR.H'
C
      IF (PSI .LT. -1.0D-14) THEN
         BETA=ABS((PSI/VGPSAT(INOD)))**VGN(INOD)
         B1=BETA+1.0D0
         B1R=1.0D0/B1
         FVGDDS=VGN1(INOD)*(BETA/PSI)*(1.0D0/PSI)*
     1               ((1.0D0 + VGN(INOD)*(BETA - 1.0D0))/
     2               (ABS(B1)**VGM(INOD)))*B1R*B1R
      ELSE
         FVGDDS=0.0D0
      END IF
C
      RETURN
      END
