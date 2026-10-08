C
C**************************  ANIS3D ************************************
C
C  rotate the current element so that the local x-axis is aligned with 
C  the velocity vector
C
C***********************************************************************
C
      SUBROUTINE ANIS3D(UUL,VVL,WWL,XX,YY,ZZ,VEL)
C
      IMPLICIT  NONE
      INTEGER   I
      REAL*8    UMOD,WMOD
      REAL*8    UUL,VVL,WWL,VEL
      REAL*8    XA(4),YA(4),ZA(4)
      REAL*8    XX(4),YY(4),ZZ(4)
C
      VEL =DSQRT(UUL*UUL + VVL*VVL + WWL*WWL)
      IF (VEL .EQ. 0.0D0) RETURN
      UMOD=DSQRT(VVL*VVL + UUL*UUL)
      WMOD=UMOD*VEL
      
      DO I=1,4
         XA(I)=XX(I)
         YA(I)=YY(I)
         ZA(I)=ZZ(I)
      END DO
      IF (UMOD .EQ. 0.0D0) THEN
         DO I=1,4
            XX(I)=ZA(I)
            YY(I)=XA(I)
            ZZ(I)=YA(I)
         END DO
      ELSE
         DO I=1,4
            XX(I)=(UUL*XA(I) + VVL*YA(I) + WWL*ZA(I))/VEL
            YY(I)=(VVL*XA(I) - UUL*YA(I))/UMOD
            ZZ(I)=(UUL*WWL*XA(I) + VVL*WWL*YA(I) -
     1            (UUL*UUL + VVL*VVL)*ZA(I))/WMOD
         END DO
      END IF
C
      RETURN
      END
