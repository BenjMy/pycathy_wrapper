C
C**************************  INIT0I ************************************
C
C  initialize an integer array to 0
C
C***********************************************************************
C
      SUBROUTINE INIT0I(NUM,IVEC)
C
      IMPLICIT  NONE
      INTEGER   I
      INTEGER   NUM
      INTEGER   IVEC(*)
C
      DO I=1,NUM
         IVEC(I)=0
      END DO
C
      RETURN
      END
