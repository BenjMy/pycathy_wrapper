C
C**************************  FHUDK2 ************************************
C
C  calculate derivative of relative hydraulic conductivity with respect
C  to pressure head using Huyakorn characteristic equation, 
C  with Kr=Se**n relationship
C
C***********************************************************************
C
      DOUBLE PRECISION FUNCTION FHUDK2(PSI,SE,DSEDP)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      REAL*8   PSI,SE,DSEDP
      INCLUDE 'SOILCHAR.H'
C
      IF (PSI .LT. HUPSIA) THEN
         FHUDK2=HUN*(SE**HUN1)*DSEDP
      ELSE
         FHUDK2=0.0D0
      END IF
C
      RETURN
      END
