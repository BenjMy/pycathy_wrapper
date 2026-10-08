C
C**************************  INIT0  ************************************
C
C  initialize a real array to 0.0
C
C***********************************************************************
C
      SUBROUTINE INIT0(N,A)
C
      IMPLICIT  NONE
      INTEGER   I
      INTEGER   N
      REAL*8    A(*)
C
      DO I=1,N
         A(I)=0.0D0
      END DO
C
      RETURN
      END
