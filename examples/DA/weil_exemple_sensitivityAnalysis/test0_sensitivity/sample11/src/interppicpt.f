C
C**************************  INTERPPICPT *******************************
C
C  use the lookup table to interpolate moisture curve values
C  (Picard and peat soil case)
C
C***********************************************************************
C
      SUBROUTINE INTERPPICPT(NLKP,ISTR,MTYPE,PAVG,KRLOC,ETALOC,SELOC)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I
      INTEGER  NLKP,ISTR,MTYPE
      REAL*8   PAVG,KRLOC,ETALOC,SELOC
      INCLUDE 'SOILCHAR.H'
C
      IF (PAVG .GT. -1.0D-14) THEN
         KRLOC=1.0D0
         SELOC=1.0D0
         ETALOC=0.0D0
      ELSE IF (PAVG .LE. PCAP(ISTR,MTYPE,1)) THEN
         KRLOC=KRWC(ISTR,MTYPE,1)
         SELOC=SATC(ISTR,MTYPE,1)
         ETALOC=(SATC(ISTR,MTYPE,2) - SELOC)/
     1          (PCAP(ISTR,MTYPE,2) - PCAP(ISTR,MTYPE,1))
      ELSE IF (PAVG .GE. PCAP(ISTR,MTYPE,NLKP)) THEN
         KRLOC=KRWC(ISTR,MTYPE,NLKP)
         SELOC=SATC(ISTR,MTYPE,NLKP)
         ETALOC=(SELOC - SATC(ISTR,MTYPE,NLKP-1))/
     1          (PCAP(ISTR,MTYPE,NLKP) - PCAP(ISTR,MTYPE,NLKP-1))
      ELSE
         I=1
         DO WHILE ((PAVG .LT. PCAP(ISTR,MTYPE,I))  .OR.
     1             (PAVG .GE. PCAP(ISTR,MTYPE,I+1)))
            I=I+1
         END DO
         KRLOC=KRWC(ISTR,MTYPE,I) +
     1         (KRWC(ISTR,MTYPE,I+1) - KRWC(ISTR,MTYPE,I))/
     2         (PCAP(ISTR,MTYPE,I+1) - PCAP(ISTR,MTYPE,I)) *
     3         (PAVG - PCAP(ISTR,MTYPE,I))
         SELOC=SATC(ISTR,MTYPE,I) +
     1         (SATC(ISTR,MTYPE,I+1) - SATC(ISTR,MTYPE,I))/
     2         (PCAP(ISTR,MTYPE,I+1) - PCAP(ISTR,MTYPE,I)) *
     3         (PAVG - PCAP(ISTR,MTYPE,I))
         ETALOC=(SELOC - SATC(ISTR,MTYPE,I))/(PAVG - PCAP(ISTR,MTYPE,I))
      END IF
C
      RETURN
      END
