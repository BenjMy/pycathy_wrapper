C
C************************** TRANSFER_SURF_F3D *****************************
C transfers the information from SURF_ROUTE to FLOW3D. The two modules are 
C based on different cells index. Then to allow the exchange of information 
C the following procedure is applied:
C
C case 1) NO "no-data" in DEM matrix
C     Surface index              Subsurface index
C     (I_BASIN_SURF)               (I_BASIN_SUB)
C  |----|----|----|---|        |----|----|----|---|
C  | 3  | 6  | 9  |12 |        | 1  | 2  | 3  |4  |
C  |----|----|----|---|        |----|----|----|---|
C  | 2  | 5  | 8  |11 |  ---\  | 5  | 6  | 7  |8  |
C  |----|----|----|---|  ---/  |----|----|----|---|
C  | 1  | 4  | 7  |10 |        | 9  | 10 | 11 |12 |
C  |----|----|----|---|        |----|----|----|---| 
C 
C case 2) "no-data" in DEM matrix
C     Surface index              Subsurface index
C     (I_BASIN_SURF)               (I_BASIN_SUB)
C  |----|----|----|---|        |----|----|----|---|
C  | 3  | ND | 9  |ND |        | 1  | ND | 3  |ND |
C  |----|----|----|---|        |----|----|----|---|
C  | 2  | 5  | 8  |11 |  ---\  | 5  | 6  | 7  |8  |
C  |----|----|----|---|  ---/  |----|----|----|---|
C  | ND | ND | 7  |10 |        | ND | ND | 11 |12 |
C  |----|----|----|---|        |----|----|----|---| 
C
C ND = No Data 
C***********************************************************************
C

      SUBROUTINE TRANSFER_SURF_F3D(NROW,NCOL,INDEX,INDEX_WITH_LAKES,
     1                           NCELL_WITH_LAKES,H_WATER_KKP1_SN,
     2                           PONDCEL,H_POOL_KKP1_VEC,
     3                           LAKES_MAP,ELEVATION_WITH_LAKES)

      IMPLICIT NONE
      INCLUDE 'CATHY.H'

      INTEGER NROW, NCOL,IROW,ICOL,I,J,JR
      INTEGER NCELL_WITH_LAKES,LAKE_NUMBER
      INTEGER INDEX(ROWMAX,*),INDEX_WITH_LAKES(ROWMAX,*)
      INTEGER LAKES_MAP(ROWMAX,*)
      REAL*8  PONDCEL(*),H_WATER_KKP1_SN(*)
      REAL*8  H_POOL_KKP1_VEC(*)
      REAL*8 ELEVATION_WITH_LAKES(*)
      INTEGER I_BASIN_SURF,I_BASIN_SUB

      
      DO ICOL=1,NCOL
         DO IROW=1,NROW
            I_BASIN_SUB=NCOL*(NROW-IROW)+ICOL
            I_BASIN_SURF=(ICOL-1)*NROW + IROW

            JR=MOD(I_BASIN_SUB,NCOL)
            IF (JR.NE.0) THEN
               J=JR
               I=(I_BASIN_SUB-J)/NCOL+1
            ELSE
               J=NCOL
               I=I_BASIN_SUB/NCOL
            END IF
            
            IF((INDEX(I,J).NE.0).AND.
     &           (INDEX_WITH_LAKES(I,J).NE.0))THEN
               
               PONDCEL(I_BASIN_SUB) = H_WATER_KKP1_SN(I_BASIN_SURF)
            END IF

            IF(LAKES_MAP(I,J).GT.0) THEN
               
               LAKE_NUMBER=LAKES_MAP(IROW,ICOL)
               PONDCEL(INDEX_WITH_LAKES(IROW,ICOL))=
     &              H_POOL_KKP1_VEC(LAKE_NUMBER)-
     &              ELEVATION_WITH_LAKES(INDEX_WITH_LAKES(IROW,ICOL))
               
               IF(PONDCEL(INDEX_WITH_LAKES(IROW,ICOL)).LT.0.D0)THEN
                  PONDCEL(INDEX_WITH_LAKES(IROW,ICOL))=0.0D0
               END IF
 
            END IF
            
         END DO

      END DO

      RETURN
      END
      
