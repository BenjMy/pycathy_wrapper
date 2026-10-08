C
C************************** TRANSFER_F3D_SURF ******************************
C transfers the information from FLOW3D to SURF_ROUTE. The two modules are 
C based on different cells index. Then to allow the exchange of information 
C the following procedure is applied:
C
C case 1) NO "no-data" in DEM matrix
C    Subsurface index             Surface index
C     (I_BASIN_SUB)               (I_BASIN_SURF)
C  |----|----|----|---|        |----|----|----|---|
C  | 1  | 2  | 3  |4  |        | 3  | 6  | 9  |12 |
C  |----|----|----|---|        |----|----|----|---|
C  | 5  | 6  | 7  |8  |  ---\  | 2  | 5  | 8  |11 |
C  |----|----|----|---|  ---/  |----|----|----|---|
C  | 9  | 10 | 11 |12 |        | 1  | 4  | 7  |10 |
C  |----|----|----|---|        |----|----|----|---|  
C case 2) "no-data" in DEM matrix
C    Subsurface index             Surface index
C     (I_BASIN_SUB)               (I_BASIN_SURF)
C  |----|----|----|---|        |----|----|----|---|
C  | 1  | ND | 3  |ND |        | 3  | ND | 9  |ND |
C  |----|----|----|---|        |----|----|----|---|
C  | 5  | 6  | 7  |8  |  ---\  | 2  | 5  | 8  |11 |
C  |----|----|----|---|  ---/  |----|----|----|---|
C  | ND | ND | 11 |12 |        | ND | ND | 7  |10 |
C  |----|----|----|---|        |----|----|----|---| 
C
C ND = No Data
C***************************************************************************
C

      SUBROUTINE TRANSFER_F3D_SURF(NROW,NCOL,INDEX,INDEX_WITH_LAKES,
     1                         SURFACE_WATER_SN,SURFACE_WATER_F3D,
     2                         NUM_R,LAKES_MAP)

      IMPLICIT NONE
      INCLUDE 'CATHY.H'
 
      INTEGER NROW, NCOL,IROW,ICOL,I,J,NUM_RR,JR
      INTEGER INDEX(ROWMAX,*),INDEX_WITH_LAKES(ROWMAX,*)
      INTEGER LAKES_MAP(ROWMAX,*)
      INTEGER NUM_R(*)
      INTEGER I_BASIN_SUB,I_BASIN_SURF
      REAL*8  SUM
      REAL*8  SURFACE_WATER_F3D(*),SURFACE_WATER_SN(*)


      DO ICOL=1,NCOL
         DO IROW=1,NROW
            I_BASIN_SUB=NCOL*(NROW-IROW)+ICOL
            I_BASIN_SURF=(ICOL-1)*NROW+IROW 

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
               
               SURFACE_WATER_SN(I_BASIN_SURF)=
     &              SURFACE_WATER_F3D(I_BASIN_SUB)

            END IF
                        
            IF(INDEX(I,J).NE.0)THEN
               IF(NUM_R(INDEX(I,J)).NE.0)THEN
                  SUM=SURFACE_WATER_F3D(I_BASIN_SUB)
                  NUM_RR=NUM_R(INDEX(I,J))
                  DO I=1,NROW
                     DO J=1,NCOL
                        IF(LAKES_MAP(I,J).EQ.NUM_RR)THEN
                        SUM=SUM+SURFACE_WATER_F3D(INDEX_WITH_LAKES(I,J))
                        END IF
                      END DO
                   END DO 
                   SURFACE_WATER_SN(I_BASIN_SURF)=SUM
                END IF
             END IF
          END DO
       END DO
C

      RETURN
      END
