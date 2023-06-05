C
C**************************  TIMUPD ************************************
C
C  update TIME and time step size for the next time step
C
C***********************************************************************
C
      SUBROUTINE TIMUPD(TIMEP,TIME,DELTAT,DTMIN,DTMAX,TMAX,
     1                  DTMAGA,DTMAGM,DTREDS,DTREDM,DTGMIN,
     2                  ITER,ITUNS1,ITUNS2)
C
      IMPLICIT  NONE
      INTEGER   ITER,ITUNS1,ITUNS2
      LOGICAL   DTGMIN
      REAL*8    TIMEP,TIME,DELTAT,DTMIN,DTMAX,TMAX
      REAL*8    DTMAGA,DTMAGM,DTREDS,DTREDM
C
      TIMEP=TIME
      IF (ITER .LT. ITUNS1) THEN
         DELTAT=DELTAT*DTMAGM + DTMAGA
         IF (DELTAT .GT. DTMAX) DELTAT=DTMAX
      END IF
      IF (ITER .GE. ITUNS2) THEN
         DELTAT=DELTAT*DTREDM - DTREDS
         IF (DELTAT .LT. DTMIN) DELTAT=DTMIN
      END IF
      IF ((TIME+DELTAT) .GE. TMAX) THEN
         DELTAT=TMAX - TIME
         TIME=TMAX
      ELSE 
         IF((TIME+2*DELTAT).GT.TMAX) DELTAT=(TMAX-TIME)/2
         TIME=TIME + DELTAT
      END IF
      IF (DELTAT .LE. DTMIN) THEN
         DTGMIN=.FALSE.
      ELSE
         DTGMIN=.TRUE.
      END IF
C
      RETURN
      END
