C
C**************************  MOISTAB ***********************************
C
C  calculate general storage term needed for mass balance calculations
C  by interpolating the moisture tables
C
C***********************************************************************
C
      SUBROUTINE MOISTAB(NLKP,NT,NTRI,TETRA,PNEW,PNODI,SNODI,ETAE)
C
      IMPLICIT NONE
      INTEGER  IEL,ISTR,IR,MTYPE,I,INOD
      INTEGER  NLKP,NT,NTRI
      INTEGER  TETRA(5,*)
      REAL*8   PAVG,SSAVG,FIAVG,ETALOC
      REAL*8   PNEW(*),PNODI(*),SNODI(*),ETAE(*)
      real*8   zero,quarter
      parameter (zero=0.0d0,quarter=0.25d0)
C
      DO IEL=1,NT
         ISTR=1+IEL/(NTRI*3)
         IR=MOD(IEL,NTRI*3)
         IF (IR .EQ. 0) ISTR=ISTR-1
         MTYPE=TETRA(5,IEL)
         PAVG = ZERO
         SSAVG = ZERO
         FIAVG = ZERO
         DO I=1,4
            INOD = TETRA(I,IEL)
            PAVG = PAVG + PNEW(INOD)
            SSAVG = SSAVG + SNODI(INOD)
            FIAVG = FIAVG + PNODI(INOD)
         END DO
         PAVG = PAVG*QUARTER
         SSAVG = SSAVG*QUARTER
         FIAVG = FIAVG*QUARTER
         CALL INTERP(NLKP,ISTR,MTYPE,PAVG,SSAVG,FIAVG,ETALOC)
         ETAE(IEL)=ETALOC
      END DO
C
      RETURN
      END
