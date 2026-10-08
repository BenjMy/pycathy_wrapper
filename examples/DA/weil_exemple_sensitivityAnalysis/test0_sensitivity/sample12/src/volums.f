C
C**************************  VOLUMS ************************************
C
C  calculate the volume of an element
C
C***********************************************************************
C
      SUBROUTINE VOLUMS(IP4,TETRAL,X,Y,Z,VOL)
C
      IMPLICIT  NONE
      INTEGER   I,J,M,II,NN,LC,LR
      INTEGER   IP4(4,4),IP3(3,3),TETRAL(5)
      REAL*8    A2,A3
      REAL*8    VOL
      REAL*8    X(*),Y(*),Z(*)
      REAL*8    AMEN(4)
      DATA      AMEN/-1.0D0,1.0D0,-1.0D0,1.0D0/
C
      VOL=0.0D0
      DO NN=1,4
         DO II=2,4
            IP3(II-1,1)=IP4(II,NN)
         END DO
         DO LC=2,3
            DO 7 LR=1,3
               IF(LR.EQ.3) IP3(LR,LC)=IP3(1,LC-1)
               IF(LR.EQ.3) GO TO 7
               IP3(LR,LC)=IP3(LR+1,LC-1)
    7       CONTINUE
         END DO
         A2=0.0D0
         A3=0.0D0
         DO II=1,3
            I=TETRAL(IP3(II,1))
            J=TETRAL(IP3(II,2))
            M=TETRAL(IP3(II,3))
            A3=Y(I)*Z(J) + A3
            A2=Y(I)*Z(M) + A2
         END DO
         VOL=VOL + X(TETRAL(NN))*AMEN(NN)*(A3-A2)/6.0D0
      END DO
C
      RETURN
      END
