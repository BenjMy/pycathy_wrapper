C
C**************************  STORCAL ***********************************
C
C  calculate volume of water in the subsurface
C
C***********************************************************************
C
      SUBROUTINE STORCAL(N,STORE1,STORE2,DSTORE,SW,PNODI,VOLNOD)
C
      IMPLICIT  NONE
      INTEGER   I
      INTEGER   N
      REAL*8    STORE1,STORE2,DSTORE
      REAL*8    SW(*),PNODI(*),VOLNOD(*)
C
      STORE1 = 0.0D0
      DO I=1,N
         STORE1 = STORE1 + SW(I)*VOLNOD(I)*PNODI(I)
      END DO
      STORE2 = STORE2 + DSTORE
C
      RETURN
      END
