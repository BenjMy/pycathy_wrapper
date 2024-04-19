C
C**************************  TOPIA  ************************************
C
C  calcola l'indice di riga IA da TOPOL
C
C***********************************************************************
C
      SUBROUTINE TOPIA(N,TOPOL,IA)
C
      IMPLICIT NONE
      INTEGER  I,J,K,M
      INTEGER  N
      INTEGER  TOPOL(*),IA(*)
C
      DO K=1,N
         I=TOPOL(K)
         J=TOPOL(K+1)-1
         DO M=I,J
            IA(M)=K
         END DO
      END DO
C
      RETURN
      END
