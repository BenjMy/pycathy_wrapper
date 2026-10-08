C
C******************************  EFFONE ********************************
C
C read parameters and arrays for first time step for the case of surface
C routing only. The subsurface-surface exchange term is directly read from
C file effraininp (IIN22).  
C
C***********************************************************************
C
      SUBROUTINE EFFONE(NCELNL,HSPEFF,HTIEFF,DELTA_X,DELTA_Y,TIME,
     1     EFFTIM,QOI_SN,SURFACE_WATER_INP,SURFACE_WATER_SN)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER I,J,N
      INTEGER NCELNL,HSPEFF,HTIEFF
      REAL*8  DELTA_X,DELTA_Y
      REAL*8  TIME,TIMEIN
      INTEGER QOI_SN(MAXCEL)
      REAL*8  EFFTIM(*)
      REAL*8  SURFACE_WATER_INP(2,MAXCEL)
      REAL*8  SURFACE_WATER_SN(MAXCEL)
      
      INCLUDE 'IOUNITS.H'

      EFFTIM(1)=0.0D0
      EFFTIM(2)=0.0D0

      DO I=1,MAXCEL
         SURFACE_WATER_INP(1,I)=0.0D0
         SURFACE_WATER_INP(2,I)=0.0D0
         SURFACE_WATER_SN(I)=0.0D0
      END DO
      
      READ(IIN22,*) HSPEFF,HTIEFF
      READ(IIN22,*) TIMEIN
      EFFTIM(2)=TIMEIN
      IF (HSPEFF .EQ. 0) THEN
         DO I=1,NCELNL
            J=QOI_SN(I)
            READ(IIN22,*)SURFACE_WATER_INP(2,J)
         END DO
      ELSE IF (HSPEFF .EQ. 1) THEN
         J=QOI_SN(1)
         READ(IIN22,*) SURFACE_WATER_INP(2,J)
         DO I=2,NCELNL
            J=QOI_SN(I)
            SURFACE_WATER_INP(2,J)=SURFACE_WATER_INP(2,QOI_SN(1))
         END DO
      END IF
   
      DO I=1,NCELNL
         J=QOI_SN(I)
         SURFACE_WATER_SN(J)=SURFACE_WATER_INP(1,J)*(DELTA_X*DELTA_Y)
      END DO
  
 800  RETURN
      END
