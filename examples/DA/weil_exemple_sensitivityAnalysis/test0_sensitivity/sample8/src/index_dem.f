C
C**************************  INDEX_DEM *********************************
C
C  builds the index matrix INDCELL that gives a cell number to all 
C  nonzero  cell values. 
C
C***********************************************************************
C
C  Note: The cell numbering procedure was modified according to 
C        the system implemented in the pre_processor. In particular,
C        compared to the previous version of the code (Anna Chiara version)
C        in counting the cells we number the no-data cells:
C
C    Present numbering           Previous numbering
C         system                      system
C  |----|----|----|---|        |----|----|----|---|
C  | 1  | ND | 3  |ND |        | 1  | ND | 2  |ND |
C  |----|----|----|---|        |----|----|----|---|
C  | 5  | 6  | 7  |8  |  ---\  | 3  | 4  | 5  | 6 |
C  |----|----|----|---|  ---/  |----|----|----|---|
C  | ND | ND | 11 |12 |        | ND | ND | 8  | 9 |
C  |----|----|----|---|        |----|----|----|---| 
C
C ND = No Data  
C MS
C**********************************************************************
      SUBROUTINE INDEX_DEM(ROWMAX,MAXCEL,NROW,NCOL,NCELL_NO_LAKES,
     1                     NCELL_WITH_LAKES,ZONE,INDCELL,
     2                     CELLCOL,CELLROW,LAKES_MAP,INDCELL_WITH_LAKES,
     3                     CELLCOL_WL,CELLROW_WL)

      IMPLICIT NONE

      INTEGER ROWMAX,MAXCEL,NROW,NCOL
      INTEGER NCELL_NO_LAKES,NCELL_WITH_LAKES
      INTEGER IROW,ICOL,N,I,K
      INTEGER INDCELL(ROWMAX,NCOL)
      INTEGER LAKES_MAP(ROWMAX,NCOL)
      INTEGER CELLCOL(MAXCEL),CELLROW(MAXCEL)
      INTEGER CELLCOL_WL(MAXCEL),CELLROW_WL(MAXCEL)
      INTEGER INDCELL_WITH_LAKES(ROWMAX,NCOL)
      INTEGER ZONE(ROWMAX,NCOL)

      INCLUDE 'IOUNITS.H'

C  first step: numbering of dem cells. In this first step
C  the lake cells are not counted in the construction of
C  INDCELL index. For the sake of convenience I read the 
C  zone raster input file

      N = 0
      K = 0
      DO IROW = 1,NROW
         DO ICOL = 1,NCOL
ccc   write(iout40,*) 'irow=',irow,'icol=',icol
            N = N + 1
            IF((ZONE(IROW,ICOL).NE.0).AND.
     &           (LAKES_MAP(IROW,ICOL).LE.0))THEN
               K = K + 1
               INDCELL(IROW,ICOL) = N
               CELLROW(K)=IROW
               CELLCOL(K)=ICOL
ccc   write(6,*) 'nonnull. irow,icol=',irow,icol,'dem=',
ccc     1              dem_map(irow,icol),'lakes=',lakes_map(irow,icol)
            ELSE
               INDCELL(IROW,ICOL) = 0
ccc   write(6,*) 'null. irow,icol=',irow,icol,'dem=',
ccc   1              dem_map(irow,icol),'lakes=',lakes_map(irow,icol)
            END IF
         END DO
      END DO
      NCELL_NO_LAKES=K

C  second step: re-reading of zone matrix considering in the numbering 
C  procedure the lake cells,i.e., construction of INDCELL_WITH_LAKES index

      N = 0
      K = 0
      DO IROW=1,NROW
         DO ICOL=1,NCOL
            N = N + 1
            IF(ZONE(IROW,ICOL).NE.0) THEN
               K = K + 1
               INDCELL_WITH_LAKES(IROW,ICOL)=N
               CELLROW_WL(K)=IROW
               CELLCOL_WL(K)=ICOL
            ELSE
               INDCELL_WITH_LAKES(IROW,ICOL)=0
            END IF
         END DO
      END DO
      NCELL_WITH_LAKES=K
 
c  writing in net.ris output file (IOUT40)

      WRITE(IOUT40,*) 'NCELL_NO_LAKES=',NCELL_NO_LAKES
      WRITE(IOUT40,*) 'NCELL_WITH_LAKES=',NCELL_WITH_LAKES
      WRITE(IOUT40,*) 'INDCELL(IROW,ICOL)'
      DO IROW=1,NROW
         WRITE(IOUT40,*) (INDCELL(IROW,ICOL),ICOL=1,NCOL)
      END DO

      WRITE(IOUT40,*) 'INDCELL_WITH_LAKES(IROW,ICOL)'
      DO IROW=1,NROW
         WRITE(IOUT40,*) (INDCELL_WITH_LAKES(IROW,ICOL),ICOL=1,NCOL)
      END DO

      WRITE(IOUT40,*) 'CELLCOL   CELLROW'
      DO I=1,NCELL_NO_LAKES
      WRITE(IOUT40,100) CELLCOL(I),CELLROW(I)
      END DO

      
100   FORMAT(I2,5X,I2)
 
      RETURN
      END
