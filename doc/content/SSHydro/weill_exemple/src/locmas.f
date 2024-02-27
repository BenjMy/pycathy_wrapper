C
C**************************  LOCMAS ************************************
C
C  set up LMASS, the part of the local mass matrix which is constant
C  for all elements (i.e. without the storage coefficient and without
C  the volume term)
C  per il caso mass lumping si procede come segue:
C  - si sommano i coefficienti extra diagonali al coefficiente diagonale
C  - si azzerano tutti i coefficienti extra diagonali
C
C***********************************************************************
C
      SUBROUTINE LOCMAS(LMASS,LUMP)
C
      IMPLICIT  NONE
      INTEGER   I,K
      INTEGER   LUMP
      REAL*8    R4,R10,R20
      REAL*8    LMASS(4,4)
C
      R4=1.0D0/4.0D0
      R10=1.0D0/10.0D0
      R20=1.0D0/20.0D0
      DO K=1,4
         DO I=K,4
            IF (LUMP .EQ. 0) THEN
               IF (K.EQ.I) THEN
                  LMASS(K,I)=R10
               ELSE
                  LMASS(K,I)=R20
               END IF
            ELSE
               IF (K.EQ.I) THEN
                  LMASS(K,I)=R4
               ELSE
                  LMASS(K,I)=0.0D0
               END IF
            END IF
         END DO
      END DO
      DO K=2,4
         DO I=1,K-1
            LMASS(K,I)=LMASS(I,K)
         END DO
      END DO
C
      RETURN
      END
