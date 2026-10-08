C
C**************************  VCOPYI ************************************
C
C  copy integer vector IVEC2 into integer vector IVEC1
C
C***********************************************************************
C
      SUBROUTINE VCOPYI(NUM,IVEC1,IVEC2)
C
      IMPLICIT  NONE
      INTEGER   I
      INTEGER   NUM
      INTEGER   IVEC1(*),IVEC2(*)
C
      DO I=1,NUM
         IVEC1(I)=IVEC2(I)
      END DO
C
      RETURN
      END
