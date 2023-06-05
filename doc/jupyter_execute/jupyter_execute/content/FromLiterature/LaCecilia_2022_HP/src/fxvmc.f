C
C**************************  FXVMC *************************************
C
C  calculate moisture content using extended van Genuchten 
C  characteristic equation
C
C***********************************************************************
C
      DOUBLE PRECISION FUNCTION FXVMC(PSI,SS,POR,INOD)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  INOD
      REAL*8   PNOT,TSR,BETA,B01
      REAL*8   PSI,SS,POR
      INCLUDE 'SOILCHAR.H'
C
      PNOT=VGPNOT(INOD)
      TSR=POR-VGRMC(INOD)
      IF (PSI .LT. PNOT) THEN
         BETA=(PSI/VGPSAT(INOD))**VGN(INOD)
         FXVMC=VGRMC(INOD) + (TSR/((BETA+1.0D0)**VGM(INOD)))
      ELSE
         B01=((PNOT/VGPSAT(INOD))**VGN(INOD)) + 1.0D0
         FXVMC=VGRMC(INOD) + TSR*(B01**(-VGM(INOD))) + SS*(PSI-PNOT)
      END IF
C
      RETURN
      END
