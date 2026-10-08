C
C**************************  AREAS  ************************************
C
C  calculate the area of an element
C
C***********************************************************************
C
      SUBROUTINE AREAS(IP3,TRIANL,X,Y,ARE)
C
      IMPLICIT  NONE
      INTEGER   I,J,M,II
      INTEGER   IP3(3,3),TRIANL(4)
      REAL*8    A2,A3
      REAL*8    ARE
      REAL*8    X(*),Y(*)
C
      A2=0.0D0
      A3=0.0D0
      DO II=1,3
         I=TRIANL(IP3(II,1))
         J=TRIANL(IP3(II,2))
         M=TRIANL(IP3(II,3))
         A3=X(I)*Y(J) + A3
         A2=X(I)*Y(M) + A2
      END DO
      ARE=0.5D0*(A3-A2)
C
      RETURN
      END
