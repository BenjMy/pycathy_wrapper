C
C**************************  WTDEPTH ***********************************
C
C  output of the water table depth at the NODVP nodes
C
C***********************************************************************
C
      SUBROUTINE WTDEPTH(NUMVP,NODVP,NSTR,NNOD,TIME,Z,PNEW)
C
      IMPLICIT  NONE
      INTEGER   I,J,INOD1,INOD2
      INTEGER   NNOD,NUMVP,NSTR
      INTEGER   NODVP(*),FLAG(NUMVP)
      REAL*8    TIME,ZERO,RC
      REAL*8    PNEW(*),Z(*),WT(NUMVP)
      INCLUDE  'IOUNITS.H'
C
      ZERO=0.0d0
      DO I=1,NUMVP
         FLAG(I)=0
CM       WT(I)=ZERO
         WT(I)=Z(NODVP(I))
      END DO
      DO I=1,NUMVP
         DO J=NSTR,1,-1
            INOD1=NODVP(I)+J*NNOD
            INOD2=NODVP(I)+(J-1)*NNOD
            IF (PNEW(INOD1).GE.ZERO.AND.PNEW(INOD2).LT.ZERO
     &         .AND.FLAG(I).EQ.0) THEN
               RC=(Z(INOD1)-Z(INOD2))/(PNEW(INOD1)-PNEW(INOD2))
CM             WT(I)=Z(NODVP(I))-(Z(INOD1)-RC*PNEW(INOD1))
               WT(I)=Z(INOD1)-RC*PNEW(INOD1)
               FLAG(I)=1
            ELSE IF (PNEW(INOD1).GE.ZERO.AND.PNEW(INOD2).LT.ZERO
     &         .AND.FLAG(I).EQ.1) THEN
               FLAG(I)=2
            ELSE IF (J.EQ.1.AND.PNEW(INOD2).GE.ZERO
     &         .AND.FLAG(I).EQ.0) THEN
               FLAG(I)=3
               WT(I)=Z(NODVP(I))+PNEW(INOD2)
            ELSE IF (J.EQ.1.AND.FLAG(I).EQ.0) THEN
               FLAG(I)=4
CM             WT(I)=Z(NODVP(I))-Z(NODVP(I)+NSTR*NNOD)
               WT(I)=Z(NODVP(I)+NSTR*NNOD)
            END IF
         END DO
      END DO
      WRITE(IOUT57,1020) TIME,(WT(I),I=1,NUMVP)
C
      RETURN
 1020 FORMAT(70(1PE15.6))
      END
