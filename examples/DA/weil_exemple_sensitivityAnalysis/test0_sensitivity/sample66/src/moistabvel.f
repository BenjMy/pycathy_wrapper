C
C**************************  MOISTABVEL ********************************
C
C  calculate characteristic curves needed for storage and velocity
C  calculations and for output by interpolating the moisture tables
C
C***********************************************************************
C
      SUBROUTINE MOISTABVEL(NLKP,NT,NTRI,TETRA,PNEW,SWE,CKRWE)
C
      IMPLICIT NONE
      INTEGER  IEL,ISTR,IR,MTYPE,I,INOD
      INTEGER  NLKP,NT,NTRI
      INTEGER  TETRA(5,*)
      REAL*8   PAVG,KRLOC,SWLOC
      REAL*8   PNEW(*),SWE(*),CKRWE(*)
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
         CALL INTERPVEL(NLKP,ISTR,MTYPE,PAVG,KRLOC,SWLOC)
         CKRWE(IEL)=KRLOC
         SWE(IEL)=SWLOC
      END DO
C
      RETURN
      END
