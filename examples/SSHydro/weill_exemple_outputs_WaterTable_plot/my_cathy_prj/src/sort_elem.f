C
C***********************************************************************
C
C     SUBROUTINE SORT_ELEM: ordina in senso crescente un vettore di N
C     interi.
C
C***********************************************************************
C
      SUBROUTINE SORT_ELEM(N,V)
C
      IMPLICIT NONE
C
      INTEGER N,V(N)
      INTEGER I,J,KK

      DO I=1,N-1
         DO J=I+1,N
            IF (V(J).LT.V(I)) THEN
               KK   = V(I)
               V(I) = V(J)
               V(J) = KK
            ENDIF
         ENDDO
      ENDDO
C
      RETURN
C
      END      
