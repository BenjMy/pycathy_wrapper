C
C**************************  ASSIGN_DEM ********************************
C
C  assign dem values to 'elevation' and multiply each cell by 'FACTOR'
C
C***********************************************************************
C
 
      SUBROUTINE ASSIGN_DEM(NROW,NCOL,DEM,INDEX,
     1                      INDEX_WITH_LAKES,LAKES_MAP,FACTOR,
     2                      ELEVATION,ELEVATION_WITH_LAKES)

      IMPLICIT NONE
      INCLUDE 'CATHY.H'

      INTEGER  NROW,NCOL
      INTEGER  IROW,ICOL,IVAR,IVAR_WITH_LAKES
      INTEGER  INDEX(ROWMAX,*)
      INTEGER  INDEX_WITH_LAKES(ROWMAX,*),LAKES_MAP(ROWMAX,*)
      REAL*8   FACTOR
      REAL*8   ELEVATION(MAXCEL),ELEVATION_WITH_LAKES(MAXCEL)
      REAL*8   DEM(ROWMAX,*)

      DO IROW = 1,NROW
       DO ICOL = 1,NCOL
        IF((DEM(IROW,ICOL).NE. 0).AND.(LAKES_MAP(IROW,ICOL).LE.0))
     1  THEN
             IVAR = INDEX(IROW,ICOL)
             ELEVATION(IVAR) = DEM(IROW,ICOL)*FACTOR
        END IF
       
        IF(DEM(IROW,ICOL).NE.0) THEN 
             IVAR_WITH_LAKES=INDEX_WITH_LAKES(IROW,ICOL)
             ELEVATION_WITH_LAKES(IVAR_WITH_LAKES)=DEM(IROW,ICOL)*FACTOR
        END IF
       END DO
      END DO

      RETURN
      END
