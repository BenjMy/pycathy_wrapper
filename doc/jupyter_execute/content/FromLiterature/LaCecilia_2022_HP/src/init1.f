C
C**************************  INIT1  ************************************
C
C  initialization of various counters, flags, and arrays
C
C***********************************************************************
C
      SUBROUTINE INIT1(NITERT,ITLIN,ITRTOT,ITER,NSTEP,
     1                 DELTAT,DTMAX,DTMIN,DTGMIN,TIMEP,TIME,
     2                 RMAX,DTSMAL,DTBIG,DTAVG,
     3                 KPRT,KSFCV,KSFCVT,KBACKT,KBACK,KLSFAI,
     4                 N,TP,SNODI,PNODI,PORE,INDE,PNEW,PTIMEP,DEF,
     5                 NP,NODDIR,CONTP,QPOLD,
     6                 HGFLAG,SFFLAG,
     7                 CPUSUB,CPUSUB_T,CPUSURF,CPUSURF_T,CPUNL,CPUVEC)
C
      IMPLICIT NONE
      INTEGER  I,J
      INTEGER  NITERT,ITLIN,ITRTOT,ITER,NSTEP
      INTEGER  KPRT,KSFCV,KSFCVT,KBACKT,KBACK,KLSFAI,N
      INTEGER  TP(*),NODDIR(*),HGFLAG(*),SFFLAG(*)
      INTEGER  NP(*),CONTP(3,*)
      LOGICAL  DTGMIN
      REAL     CPUSUB,CPUSUB_T,CPUSURF,CPUSURF_T,CPUNL,CPUVEC(*)
      REAL*8   DELTAT,DTMAX,DTMIN,TIMEP,TIME
      REAL*8   RMAX,DTSMAL,DTBIG,DTAVG
      REAL*8   SNODI(*),PNODI(*),PORE(*),INDE(*),PNEW(*),PTIMEP(*)
      REAL*8   QPOLD(*),DEF(*)
C
      DO I=1,11
         CPUVEC(I)=0.0
      END DO
      CPUNL=0.0
      CPUSUB=0.0
      CPUSUB_T=0.0
      CPUSURF=0.0
      CPUSURF_T=0.0
C
      NITERT=0
      ITLIN=0
      ITRTOT=0
      ITER=1
      NSTEP=1
      IF (DELTAT .GT. DTMAX) DELTAT=DTMAX
      IF (DELTAT .LE. DTMIN) THEN
         DELTAT=DTMIN
         DTGMIN=.FALSE.
      ELSE
         DTGMIN=.TRUE.
      END IF
      TIMEP=0.0D0
      IF (DELTAT .GE. 1.0E+15) THEN
         TIME=0.0D0
      ELSE
         TIME=DELTAT
      END IF 
      DTSMAL=RMAX
      DTBIG=0.0D0
      DTAVG=0.0D0
      KPRT=1
      KSFCV=0
      KSFCVT=0
      KBACKT=0
      KBACK=0
      KLSFAI=0
      CALL INIT0I(N,TP)
      CALL INIT0R(N,SNODI)
      CALL INIT0R(N,PNODI)
      CALL INIT0R(N,INDE)
      CALL INIT0R(N,DEF)
      CALL INIT0R(N,PORE)
      CALL VCOPYR(N,PNEW,PTIMEP)
C     CALL VCOPYI(NP,NODDIR,CONTP)
C     CALL INIT0R(NP,QPOLD)
      DO I=1,9
         HGFLAG(I)=0
      END DO
      DO I=1,5
         SFFLAG(I)=0
      END DO
C
      RETURN
      END
