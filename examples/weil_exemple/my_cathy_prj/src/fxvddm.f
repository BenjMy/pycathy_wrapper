C
C**************************  FXVDDM ************************************
C
C  calculate second derivative of moisture content with respect to
C  pressure head using extended van Genuchten characteristic equation.
C  for the extended van Genuchten equation this gives the derivative of
C  the overall storage coefficient with respect to pressure head.
C
C***********************************************************************
C
      DOUBLE PRECISION FUNCTION FXVDDM(PSI,POR,INOD)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  INOD
      REAL*8   BETA,B1,B1R
      REAL*8   PSI,POR
      INCLUDE 'SOILCHAR.H'
C
      IF (PSI .LT. VGPNOT(INOD)) THEN
         BETA=(PSI/VGPSAT(INOD))**VGN(INOD)
         B1=BETA+1.0D0
         B1R=1.0D0/B1
         FXVDDM=VGN1(INOD)*(POR-VGRMC(INOD))*(BETA/PSI)*(1.0D0/PSI)*
     1              ((1.0D0 + VGN(INOD)*(BETA - 1.0D0))/
     2              (B1**VGM(INOD)))*B1R*B1R
      ELSE
         FXVDDM=0.0D0
      END IF
C
      RETURN
      END
