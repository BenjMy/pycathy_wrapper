C
C**************************  DTSTAT ************************************
C
C  determine largest, smallest, and average time step sizes
C
C***********************************************************************
C
      SUBROUTINE DTSTAT(TIME,DELTAT,DTBIG,TBIG,DTSMAL,TSMAL,DTAVG)
C
      IMPLICIT  NONE
      REAL*8    TIME,DELTAT,DTBIG,TBIG,DTSMAL,TSMAL,DTAVG
C
      IF (DELTAT .GT. DTBIG) THEN
         DTBIG=DELTAT
         TBIG=TIME
      END IF
      IF (DELTAT .LT. DTSMAL) THEN
         DTSMAL=DELTAT
         TSMAL=TIME
      END IF
      DTAVG=DTAVG + DELTAT
C
      RETURN
      END
