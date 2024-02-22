C
C**************************  CHTANP ************************************
C
C  calculate localized tangent slope approximations of moisture curve
C  derivatives needed for Picard scheme for the case KSLOPE=4.
C  Note that DSETAN contains tangent slope values at each node only for 
C  the case IVGHU=1; for the other IVGHU cases the tangent slope is 
C  constant for all nodes and we use DSETAN(1).
C
C***********************************************************************
C
      SUBROUTINE CHTANP(N,SNODI,PNODI,IVGHU)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I
      INTEGER  N,IVGHU
      REAL*8   FVGSE,FXVMC,FHUSE,FBCSE
      REAL*8   SEL,SER,DIFF
      REAL*8   SNODI(*),PNODI(*)
      INCLUDE 'SOILCHAR.H'
C
      IF (IVGHU .EQ. 0) THEN
         I=1
         SEL=FVGSE(PSEL,I)
         SER=FVGSE(PSER,I)
         DSETAN(1)=(SER - SEL)/(PSER - PSEL)
      ELSE IF (IVGHU .EQ. 1) THEN
         DIFF=1.0D0/(PSER - PSEL)
         DO I=1,N
            SEL=FXVMC(PSEL,SNODI(I),PNODI(I),I)
            SER=FXVMC(PSER,SNODI(I),PNODI(I),I)
            DSETAN(I)=(SER - SEL)*DIFF
         END DO
      ELSE IF (IVGHU .EQ. 2) THEN
         SEL=FHUSE(PSEL)
         SER=FHUSE(PSER)
         DSETAN(1)=(SER - SEL)/(PSER - PSEL)
      ELSE IF (IVGHU .EQ. 3) THEN
         SEL=FHUSE(PSEL)
         SER=FHUSE(PSER)
         DSETAN(1)=(SER - SEL)/(PSER - PSEL)
      ELSE
         SEL=FBCSE(PSEL)
         SER=FBCSE(PSER)
         DSETAN(1)=(SER - SEL)/(PSER - PSEL)
      END IF
C
      RETURN
      END
