C
C**************************  VNOD3D ************************************
C
C  average Darcy velocitities on each tetrahedra to obtain the Darcy
C  velocities at each node
C
C***********************************************************************
C
      SUBROUTINE VNOD3D(N,NT,TP,TETRA,UU,VV,WW,UNOD,VNOD,WNOD)
C
      IMPLICIT  NONE
      INTEGER   I,K,II,INOD
      INTEGER   N,NT
      INTEGER   TP(*),TETRA(5,*)
      REAL*8    UU(*),VV(*),WW(*),UNOD(*),VNOD(*),WNOD(*)
C
      CALL INIT0R(N,UNOD)
      CALL INIT0R(N,VNOD)
      CALL INIT0R(N,WNOD)
      DO K=1,NT
         DO II=1,4
            INOD=TETRA(II,K)
            UNOD(INOD)=UNOD(INOD) + UU(K)
            VNOD(INOD)=VNOD(INOD) + VV(K)
            WNOD(INOD)=WNOD(INOD) + WW(K)
         END DO
      END DO
      DO I=1,N
         UNOD(I)=UNOD(I)/TP(I)
         VNOD(I)=VNOD(I)/TP(I)
         WNOD(I)=WNOD(I)/TP(I)
      END DO
C
      RETURN
      END
