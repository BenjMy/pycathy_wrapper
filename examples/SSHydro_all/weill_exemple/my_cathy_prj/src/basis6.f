C
C**************************  BASIS6 ************************************
C
C  calculate the local basis function coefficients. 
C  Basis function coefficients are divided by 6.
C
C***********************************************************************
C
      SUBROUTINE BASIS6(IP4,TETRAL,X,Y,Z,AIL,BIL,CIL,DIL)
C
      IMPLICIT  NONE
      INTEGER   I,J,M,II,NN,LC,LR
      INTEGER   IP4(4,4),IP3(3,3),TETRAL(5)
      REAL*8    A2,A3
      REAL*8    X(*),Y(*),Z(*),AIL(4),BIL(4),CIL(4),DIL(4)
      REAL*8    AMEN(5)
      DATA      AMEN/-1.0D0,1.0D0,-1.0D0,1.0D0,-1.0D0/
C
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
         A3=0.0D0
         DO II=1,3
            I=TETRAL(IP3(II,1))
            J=TETRAL(IP3(II,2))
            M=TETRAL(IP3(II,3))
            A3=A3 + X(M) * (Y(I)*Z(J) - Z(I)*Y(J))
         END DO
         AIL(NN)=AMEN(NN+1)*A3/6.0D0
         A2=0.0D0
         A3=0.0D0
         DO II=1,3
            I=TETRAL(IP3(II,1))
            J=TETRAL(IP3(II,2))
            M=TETRAL(IP3(II,3))
            A3=Y(I)*Z(J) + A3
            A2=Y(I)*Z(M) + A2
         END DO
         BIL(NN)=AMEN(NN)*(A3-A2)/6.0D0
         A2=0.0D0
         A3=0.0D0
         DO II=1,3
            I=TETRAL(IP3(II,1))
            J=TETRAL(IP3(II,2))
            M=TETRAL(IP3(II,3))
            A3=X(I)*Z(J) + A3
            A2=X(I)*Z(M) + A2
         END DO
         CIL(NN)=AMEN(NN+1)*(A3-A2)/6.0D0
         A2=0.0D0
         A3=0.0D0
         DO II=1,3
            I=TETRAL(IP3(II,1))
            J=TETRAL(IP3(II,2))
            M=TETRAL(IP3(II,3))
            A3=X(I)*Y(J) + A3
            A2=X(I)*Y(M) + A2
         END DO
         DIL(NN)=AMEN(NN)*(A3-A2)/6.0D0
      END DO
C
      RETURN
      END
