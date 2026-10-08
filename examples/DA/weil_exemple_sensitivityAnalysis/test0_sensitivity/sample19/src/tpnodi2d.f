C
C**************************  TPNODI2D **********************************
C
C  ripartisco ai nodi le quantita' che sono date per 
C  elementi e che devono andare ai nodi
C
C***********************************************************************
C
      SUBROUTINE TPNODI2d(nnod,ntri,TRIANG,z,TP2d,eltria)

      IMPLICIT NONE
      include 'CATHY.H'
      INTEGER  K,II,INOD,ISTOP
      INTEGER  nnod,ntri
      INTEGER  TRIANG(4,*),TP2d(*)
      REAL*8   z(*),eltria(*)
!CMS  
!CMS      REAL*8   ZTEMP(2856)
!CMS      INTEGER  jj,kk,kz
!CMS

      do k=1,nnod
        tp2d(k)=0
        z(k)=0.0d0
      end do

      DO K=1,NTRI
         DO II=1,3
            INOD=TRIANG(II,K)
            z(INOD)=z(INOD) + eltria(K)
            TP2d(INOD)=TP2d(INOD)+1
         END DO
      END DO
      ISTOP=0
      DO K=1,NNOD
         IF (TP2d(K) .EQ. 0) THEN
            WRITE(6,1000) K
            ISTOP=1
         ELSE
            z(K)=z(K)/TP2d(K)
         END IF
      END DO
!CMS      
!CMS      OPEN(98,file='./grid2d')
!CMS      DO kk=1,NNOD
!CMS          READ(98,*) ZTEMP(kk)
!CMS          z(kk) = ZTEMP(kk)
!CMS      END DO
!CMS      CLOSE(98)
!CMS      
      IF (ISTOP .NE. 0) THEN
         STOP
      END IF
C
      RETURN
 1000 FORMAT(/,' ERROR IN MESH CONSTRUCTION:    NODE ',I7,
     1         ' IS NOT CONNECTED TO ANY ELEMENTS')
      END
