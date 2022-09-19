C
C**************************  RAST_INPUT_REAL*********************************
C
C  reads a real raster map in GRASS (ASCII) format from unit IIN
C  according to cartesian system. This subroutine is used to read real 
C  surface module input files 
C
C****************************************************************************
C
      SUBROUTINE RAST_INPUT_REAL(IIN,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1     RAST)


      IMPLICIT NONE
      INCLUDE 'CATHY.H'
     
      INTEGER  IIN,NROW,NCOL
      INTEGER  IROW,ICOL
      REAL*8  RAST(COLMAX,*)
      REAL*8 NORTH,SOUTH,EAST,WEST
      CHARACTER RDWR*80,LINE*80

c
c  reads the first four lines (north,south,east,west)
c
      READ(IIN,'(A80)') LINE
      WRITE(RDWR,'(A74)') LINE(8:80)
      READ(RDWR,*) NORTH
      READ(IIN,'(A80)') LINE
      WRITE(RDWR,'(A74)') LINE(8:80)
      READ(RDWR,*) SOUTH 
      READ(IIN,'(A80)') LINE
      WRITE(RDWR,'(A74)') LINE(8:80)
      READ(RDWR,*) EAST 
      READ(IIN,'(A80)') LINE
      WRITE(RDWR,'(A74)') LINE(8:80)
      READ(RDWR,*) WEST 
c
c  reads in NROW and NCOL
c
      READ(IIN,'(A80)') LINE
      WRITE(RDWR,'(A75)') LINE(7:80)
      READ(RDWR,*) NROW
      READ(IIN,'(A80)') LINE
      WRITE(RDWR,'(A75)') LINE(7:80)
      READ(RDWR,*) NCOL
c  reads the raster map      
      DO IROW=NROW,1,-1
         READ(IIN,*) (RAST(ICOL,IROW),ICOL=1,NCOL)
      END DO

      RETURN
      END
