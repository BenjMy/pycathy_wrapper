C
C**************************  TETNEW ************************************
C
C  set up TETJA for nonsymmetric case
C
C***********************************************************************
C
      SUBROUTINE TETNEW(NT,TETRA,JA,TOPOL,TETJA)
C
      IMPLICIT NONE
      INTEGER  I,J,II,JJ,IEL,IND
      INTEGER  NT
      INTEGER  TETRA(5,*),TETJA(4,4,*),JA(*),TOPOL(*)
C
      DO IEL=1,NT
         DO I=1,4
            II=TETRA(I,IEL)
            DO J=1,4
               JJ=TETRA(J,IEL)
               IND=TOPOL(II)-1
  100          IND=IND+1
               IF (JA(IND) .NE. JJ) GO TO 100
               TETJA(I,J,IEL)=IND
            END DO
         END DO
      END DO
C
      RETURN
      END
