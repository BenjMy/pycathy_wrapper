C
C**************************  MOISTABPT *********************************
C
C  calculate general storage term needed for mass balance calculations
C  by interpolating the moisture tables (peat soil case)
C
C***********************************************************************
C
      SUBROUTINE MOISTABPT(NLKP,NT,NTRI,TETRA,PNEW,SEELT,ETAE)
C
      IMPLICIT NONE
      INTEGER  IEL,ISTR,IR,MTYPE,I,INOD
      INTEGER  NLKP,NT,NTRI
      INTEGER  TETRA(5,*)
      REAL*8   PAVG,ETALOC,SELOC
      REAL*8   PNEW(*),SEELT(*),ETAE(*)
      real*8   zero,quarter
      parameter (zero=0.0d0,quarter=0.25d0)
C
      DO IEL=1,NT
         ISTR=1+IEL/(NTRI*3)
         IR=MOD(IEL,NTRI*3)
         IF (IR .EQ. 0) ISTR=ISTR-1
         MTYPE=TETRA(5,IEL)
         PAVG = ZERO
         DO I=1,4
            INOD = TETRA(I,IEL)
            PAVG = PAVG + PNEW(INOD)
         END DO
         PAVG = PAVG*QUARTER
         CALL INTERPPT(NLKP,ISTR,MTYPE,PAVG,ETALOC,SELOC)
         ETAE(IEL)=ETALOC
         SEELT(IEL)=SELOC
      END DO
C
      RETURN
      END
