C
C**************************  VCOPYR ************************************
C
C  copy real vector RVEC2 into real vector RVEC1
C
C***********************************************************************
C
      SUBROUTINE VCOPYR(NUM,RVEC1,RVEC2)
C
      IMPLICIT  NONE
      INTEGER   I
      INTEGER   NUM
      REAL*8    RVEC1(*),RVEC2(*)
C
      DO I=1,NUM
         RVEC1(I)=RVEC2(I)
      END DO
C
      RETURN
      END
