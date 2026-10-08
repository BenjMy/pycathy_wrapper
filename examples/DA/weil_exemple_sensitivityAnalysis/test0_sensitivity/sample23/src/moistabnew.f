C
C**************************  MOISTABNEW ********************************
C
C  calculate characteristic curves for unsaturated zone
C  by interpolating the moisture tables (Newton case)
C
C***********************************************************************
C
      SUBROUTINE MOISTABNEW(NLKP,NT,NTRI,TETRA,PTNEW,PNODI,SNODI,
     1                      SWE,CKRWE,ETAE,DCKRWE,DETAIE)
C
      IMPLICIT NONE
      INTEGER  IEL,ISTR,IR,MTYPE,I,INOD
      INTEGER  NLKP,NT,NTRI
      INTEGER  TETRA(5,*)
      REAL*8   PAVG,SSAVG,FIAVG,KRLOC,ETALOC,SWLOC,DKRLOC,DETALOC
      REAL*8   PTNEW(*),PNODI(*),SNODI(*)
      REAL*8   SWE(*),CKRWE(*),ETAE(*),DCKRWE(*),DETAIE(*)
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
            PAVG = PAVG + PTNEW(INOD)
            SSAVG = SSAVG + SNODI(INOD)
            FIAVG = FIAVG + PNODI(INOD)
         END DO
         PAVG = PAVG*QUARTER
         SSAVG = SSAVG*QUARTER
         FIAVG = FIAVG*QUARTER
         CALL INTERPNEW(NLKP,ISTR,MTYPE,PAVG,SSAVG,FIAVG,
     1                  KRLOC,ETALOC,SWLOC,DKRLOC,DETALOC)
         CKRWE(IEL)=KRLOC
         ETAE(IEL)=ETALOC
         SWE(IEL)=SWLOC
         DCKRWE(IEL)=DKRLOC
         DETAIE(IEL)=DETALOC
      END DO
C
      RETURN
      END
