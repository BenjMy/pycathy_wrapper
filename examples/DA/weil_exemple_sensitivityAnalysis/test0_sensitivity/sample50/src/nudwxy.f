C
C**************************  NUDWXY ************************************
C
C  calculate horizontal weighting function W(x,y) for nudging
C
C***********************************************************************
C
      DOUBLE PRECISION FUNCTION NUDWXY(WFLAG,XX,YY,XXO,YYO,RXY)
C
      IMPLICIT NONE
      INTEGER  WFLAG
      REAL*8   R2,D2
      REAL*8   XX,YY,XXO,YYO,RXY
C
      R2=RXY*RXY
      D2=(XX-XXO)*(XX-XXO) + (YY-YYO)*(YY-YYO)
      IF (WFLAG .EQ. 0) THEN
         IF (D2 .LE. R2) THEN
            NUDWXY=(R2-D2)/(R2+D2)
         ELSE
            NUDWXY=0.0D0
         END IF
      ELSE
         NUDWXY=DEXP(-1.0d0*(D2/R2))
      END IF
C
      RETURN
      END
