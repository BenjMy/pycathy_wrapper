C
C**************************  INTERP ************************************
C
C  use the lookup table to interpolate general storage value
C
C***********************************************************************
C
      SUBROUTINE INTERP(NLKP,ISTR,MTYPE,PAVG,SSAVG,FIAVG,ETALOC)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I
      INTEGER  NLKP,ISTR,MTYPE
      REAL*8   SWLOC
      REAL*8   PAVG,SSAVG,FIAVG,ETALOC
      INCLUDE 'SOILCHAR.H'
C
      IF (PAVG .GT. -1.0D-14) THEN
         ETALOC=SSAVG
      ELSE IF (PAVG .LE. PCAP(ISTR,MTYPE,1)) THEN
         SWLOC=SATC(ISTR,MTYPE,1)
         ETALOC=SSAVG*SWLOC +
     1          FIAVG*(SATC(ISTR,MTYPE,2) - SWLOC)/
     2                (PCAP(ISTR,MTYPE,2) - PCAP(ISTR,MTYPE,1))
      ELSE IF (PAVG .GE. PCAP(ISTR,MTYPE,NLKP)) THEN
         SWLOC=SATC(ISTR,MTYPE,NLKP)
         ETALOC=SSAVG*SWLOC +
     1          FIAVG*(SWLOC - SATC(ISTR,MTYPE,NLKP-1))/
     2                (PCAP(ISTR,MTYPE,NLKP) - PCAP(ISTR,MTYPE,NLKP-1))
      ELSE
         I=1
         DO WHILE ((PAVG .LT. PCAP(ISTR,MTYPE,I))  .OR.
     1             (PAVG .GE. PCAP(ISTR,MTYPE,I+1)))
            I=I+1
         END DO
         SWLOC=SATC(ISTR,MTYPE,I) +
     1         (SATC(ISTR,MTYPE,I+1) - SATC(ISTR,MTYPE,I))/
     2         (PCAP(ISTR,MTYPE,I+1) - PCAP(ISTR,MTYPE,I)) *
     3         (PAVG - PCAP(ISTR,MTYPE,I))
         ETALOC=SSAVG*SWLOC +
     1          FIAVG*(SWLOC - SATC(ISTR,MTYPE,I))/
     2                (PAVG - PCAP(ISTR,MTYPE,I))
      END IF
C
      RETURN
      END
