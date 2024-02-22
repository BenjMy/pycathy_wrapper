C
C**************************  ELTNOD ************************************
C
C  average generic vector on each tetrahedra to obtain the corresponding  
C  vector at each node
C
C***********************************************************************
C
      SUBROUTINE ELTNOD(N,NT,TP,TETRA,VELT,VNOD)
C
      IMPLICIT  NONE
      INTEGER   I,K,II,INOD
      INTEGER   N,NT
      INTEGER   TP(*),TETRA(5,*)
      REAL*8    VELT(*),VNOD(*)
C
      CALL INIT0R(N,VNOD)
      DO K=1,NT
         DO II=1,4
            INOD=TETRA(II,K)
            VNOD(INOD)=VNOD(INOD) + VELT(K)
         END DO
      END DO
      DO I=1,N
         VNOD(I)=VNOD(I)/TP(I)
      END DO
C
      RETURN
      END
