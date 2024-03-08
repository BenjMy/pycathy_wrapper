C
C**************************  NUDWZ  ************************************
C
C  calculate vertical weighting function W(z) for nudging
C
C***********************************************************************
C
      DOUBLE PRECISION FUNCTION NUDWZ(WFLAG,ZZ,ZZO,RZ)
C
      IMPLICIT NONE
      INTEGER  WFLAG
      REAL*8   ZD
      REAL*8   ZZ,ZZO,RZ
C
      ZD=DABS(ZZ-ZZO)
      IF (WFLAG .EQ. 0) THEN
         IF (ZD .LE. RZ) THEN
            NUDWZ=1.0D0 - ZD/RZ
         ELSE
            NUDWZ=0.0D0
         END IF
      ELSE
         NUDWZ=DEXP(-1.0d0*(ZD/RZ)**2.0d0)
      END IF
C
      RETURN
      END
