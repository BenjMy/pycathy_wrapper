C
C**************************  INIT0R ************************************
C
C  initialize a real array to 0.0
C
C***********************************************************************
C
      SUBROUTINE INIT0R(NUM,RVEC)
C
      IMPLICIT  NONE
      INTEGER   I
      INTEGER   NUM
      REAL*8    RVEC(*)
C
      DO I=1,NUM
         RVEC(I)=0.0D0
      END DO
C
      RETURN
      END
