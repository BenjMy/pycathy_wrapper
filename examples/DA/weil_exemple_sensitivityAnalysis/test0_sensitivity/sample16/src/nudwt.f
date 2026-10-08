C
C**************************  NUDWT  ************************************
C
C  calculate temporal weighting function W(t) for nudging
C
C***********************************************************************
C
      DOUBLE PRECISION FUNCTION NUDWT(WFLAG,TT,TTO,TAU)
C
      IMPLICIT NONE
      INTEGER  WFLAG
      REAL*8   TD,TAU2
      REAL*8   TT,TTO,TAU
C
      IF (WFLAG .EQ. 0) THEN
         TD=DABS(TT-TTO)
         TAU2=0.5D0*TAU
         IF (TD .LT. TAU2) THEN 
            NUDWT=1.0D0
         ELSE IF (TD .GE. TAU2  .AND.  TD .LE. TAU) THEN
            NUDWT=(TAU - TD)/TAU2
         ELSE
            NUDWT=0.0D0
         END IF
      ELSE
         TD=TT-TTO
         IF (TD .LT. 0.0D0) THEN
            NUDWT=0.0D0
         ELSE
            NUDWT=DEXP(-TD/TAU)
         END IF
      END IF
C
      RETURN
      END
