C
C**************************  MOISTABPICPT ******************************
C
C  calculate characteristic curves for unsaturated zone
C  by interpolating the moisture tables (Picard and peat soil case) 
C
C***********************************************************************
C
      SUBROUTINE MOISTABPICPT(NLKP,NT,NTRI,TETRA,PTNEW,SEELT,CKRWE,ETAE)
C
      IMPLICIT NONE
      INTEGER  IEL,ISTR,IR,MTYPE,I,INOD
      INTEGER  NLKP,NT,NTRI
      INTEGER  TETRA(5,*)
      REAL*8   PAVG,KRLOC,ETALOC,SELOC
      REAL*8   PTNEW(*),SEELT(*),CKRWE(*),ETAE(*)
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
            PAVG = PAVG + PTNEW(INOD)
         END DO
         PAVG = PAVG*QUARTER
         CALL INTERPPICPT(NLKP,ISTR,MTYPE,PAVG,KRLOC,ETALOC,SELOC)
         CKRWE(IEL)=KRLOC
         ETAE(IEL)=ETALOC
         SEELT(IEL)=SELOC
      END DO
C
      RETURN
      END
