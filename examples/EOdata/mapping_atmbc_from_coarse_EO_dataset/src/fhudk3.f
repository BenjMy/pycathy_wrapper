C
C**************************  FHUDK3 ************************************
C
C  calculate derivative of relative hydraulic conductivity with respect
C  to pressure head using Huyakorn characteristic equation, 
C  with conductivity relationship from 
C  Table 3 of 1984 paper (log_10 Kr(Se) curve)
C
C***********************************************************************
C
      DOUBLE PRECISION FUNCTION FHUDK3(PSI,SE,DSEDP,KRW)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      REAL*8   PSI,SE,DSEDP,KRW
      INCLUDE 'SOILCHAR.H'
C
      IF (PSI .LT. HUPSIA) THEN
         FHUDK3=(HU2A*SE + HUB2A)*DSEDP*KRW*HULN10
      ELSE
         FHUDK3=0.0D0
      END IF
C
      RETURN
      END
