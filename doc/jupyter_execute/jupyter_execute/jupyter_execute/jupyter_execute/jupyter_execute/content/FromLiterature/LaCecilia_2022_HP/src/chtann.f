C
C**************************  CHTANN ************************************
C
C  calculate localized tangent slope approximations of moisture curve
C  derivatives needed for Newton scheme for the case KSLOPE=4.
C  Note that DSETAN, DDSE1T, and DDSE2T contain tangent slope
C  values at each node only for the case IVGHU=1; for the other
C  IVGHU cases the tangent slope values are constant for all nodes
C  and we use DSETAN(1), DDSE1T(1), and DDSE2T(1).
C
C***********************************************************************
C
      SUBROUTINE CHTANN(N,SNODI,PNODI,IVGHU)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I
      INTEGER  N,IVGHU
      REAL*8   FVGKR,FVGSE,FVGDSE,FXVKR,FXVMC,FXVDMC
      REAL*8   FHUKR2,FHUKR3,FHUSE,FHUDSE,FBCKR,FBCSE,FBCDSE
      REAL*8   KRL,KRR,SEL,SER,DSE1L,DSE1R,DSE2L,DSE2R
      REAL*8   DIFF,DIFF1,DIFF2
      REAL*8   SNODI(*),PNODI(*)
      INCLUDE 'SOILCHAR.H'
C
      IF (IVGHU .EQ. 0) THEN
         I=1
         SEL=FVGSE(PSEL,I)
         SER=FVGSE(PSER,I)
         DSETAN(1)=(SER - SEL)/(PSER - PSEL)
         IF (PDSE1L .GE. PSEL  .AND.  PDSE1L .LE. PSER) THEN
            DSE1L=DSETAN(1)
         ELSE
            DSE1L=FVGDSE(PDSE1L,I)
         END IF
         IF (PDSE1R .GE. PSEL  .AND.  PDSE1R .LE. PSER) THEN
            DSE1R=DSETAN(1)
         ELSE
            DSE1R=FVGDSE(PDSE1R,I)
         END IF
         DDSE1T(1)=(DSE1R - DSE1L)/(PDSE1R - PDSE1L)
         IF (PDSE2L .GE. PSEL  .AND.  PDSE2L .LE. PSER) THEN
            DSE2L=DSETAN(1)
         ELSE
            DSE2L=FVGDSE(PDSE2L,I)
         END IF
         IF (PDSE2R .GE. PSEL  .AND.  PDSE2R .LE. PSER) THEN
            DSE2R=DSETAN(1)
         ELSE
            DSE2R=FVGDSE(PDSE2R,I)
         END IF
         DDSE2T(1)=(DSE2R - DSE2L)/(PDSE2R - PDSE2L)
         SEL=FVGSE(PKRL,I)
         SER=FVGSE(PKRR,I)
         KRL=FVGKR(PKRL,SEL,I)
         KRR=FVGKR(PKRR,SER,I)
         DKRTAN=(KRR - KRL)/(PKRR - PKRL)
      ELSE IF (IVGHU .EQ. 1) THEN
         DIFF=1.0D0/(PSER - PSEL)
         DIFF1=1.0D0/(PDSE1R - PDSE1L)
         DIFF2=1.0D0/(PDSE2R - PDSE2L)
         DO I=1,N
            SEL=FXVMC(PSEL,SNODI(I),PNODI(I),I)
            SER=FXVMC(PSER,SNODI(I),PNODI(I),I)
            DSETAN(I)=(SER - SEL)*DIFF
            IF (PDSE1L .GE. PSEL  .AND.  PDSE1L .LE. PSER) THEN
               DSE1L=DSETAN(I)
            ELSE
               DSE1L=FXVDMC(PDSE1L,SNODI(I),PNODI(I),I)
            END IF
            IF (PDSE1R .GE. PSEL  .AND.  PDSE1R .LE. PSER) THEN
               DSE1R=DSETAN(I)
            ELSE
               DSE1R=FXVDMC(PDSE1R,SNODI(I),PNODI(I),I)
            END IF
            DDSE1T(I)=(DSE1R - DSE1L)*DIFF1
            IF (PDSE2L .GE. PSEL  .AND.  PDSE2L .LE. PSER) THEN
               DSE2L=DSETAN(I)
            ELSE
               DSE2L=FXVDMC(PDSE2L,SNODI(I),PNODI(I),I)
            END IF
            IF (PDSE2R .GE. PSEL  .AND.  PDSE2R .LE. PSER) THEN
               DSE2R=DSETAN(I)
            ELSE
               DSE2R=FXVDMC(PDSE2R,SNODI(I),PNODI(I),I)
            END IF
            DDSE2T(I)=(DSE2R - DSE2L)*DIFF2
         END DO
         KRL=FXVKR(PKRL,I)
         KRR=FXVKR(PKRR,I)
         DKRTAN=(KRR - KRL)/(PKRR - PKRL)
      ELSE IF (IVGHU .EQ. 2) THEN
         SEL=FHUSE(PSEL)
         SER=FHUSE(PSER)
         DSETAN(1)=(SER - SEL)/(PSER - PSEL)
         IF (PDSE1L .GE. PSEL  .AND.  PDSE1L .LE. PSER) THEN
            DSE1L=DSETAN(1)
         ELSE
            DSE1L=FHUDSE(PDSE1L)
         END IF
         IF (PDSE1R .GE. PSEL  .AND.  PDSE1R .LE. PSER) THEN
            DSE1R=DSETAN(1)
         ELSE
            DSE1R=FHUDSE(PDSE1R)
         END IF
         DDSE1T(1)=(DSE1R - DSE1L)/(PDSE1R - PDSE1L)
         IF (PDSE2L .GE. PSEL  .AND.  PDSE2L .LE. PSER) THEN
            DSE2L=DSETAN(1)
         ELSE
            DSE2L=FHUDSE(PDSE2L)
         END IF
         IF (PDSE2R .GE. PSEL  .AND.  PDSE2R .LE. PSER) THEN
            DSE2R=DSETAN(1)
         ELSE
            DSE2R=FHUDSE(PDSE2R)
         END IF
         DDSE2T(1)=(DSE2R - DSE2L)/(PDSE2R - PDSE2L)
         SEL=FHUSE(PKRL)
         SER=FHUSE(PKRR)
         KRL=FHUKR2(PKRL,SEL)
         KRR=FHUKR2(PKRR,SER)
         DKRTAN=(KRR - KRL)/(PKRR - PKRL)
      ELSE IF (IVGHU .EQ. 3) THEN
         SEL=FHUSE(PSEL)
         SER=FHUSE(PSER)
         DSETAN(1)=(SER - SEL)/(PSER - PSEL)
         IF (PDSE1L .GE. PSEL  .AND.  PDSE1L .LE. PSER) THEN
            DSE1L=DSETAN(1)
         ELSE
            DSE1L=FHUDSE(PDSE1L)
         END IF
         IF (PDSE1R .GE. PSEL  .AND.  PDSE1R .LE. PSER) THEN
            DSE1R=DSETAN(1)
         ELSE
            DSE1R=FHUDSE(PDSE1R)
         END IF
         DDSE1T(1)=(DSE1R - DSE1L)/(PDSE1R - PDSE1L)
         IF (PDSE2L .GE. PSEL  .AND.  PDSE2L .LE. PSER) THEN
            DSE2L=DSETAN(1)
         ELSE
            DSE2L=FHUDSE(PDSE2L)
         END IF
         IF (PDSE2R .GE. PSEL  .AND.  PDSE2R .LE. PSER) THEN
            DSE2R=DSETAN(1)
         ELSE
            DSE2R=FHUDSE(PDSE2R)
         END IF
         DDSE2T(1)=(DSE2R - DSE2L)/(PDSE2R - PDSE2L)
         SEL=FHUSE(PKRL)
         SER=FHUSE(PKRR)
         KRL=FHUKR3(PKRL,SEL)
         KRR=FHUKR3(PKRR,SER)
         DKRTAN=(KRR - KRL)/(PKRR - PKRL)
      ELSE
         SEL=FBCSE(PSEL)
         SER=FBCSE(PSER)
         DSETAN(1)=(SER - SEL)/(PSER - PSEL)
         IF (PDSE1L .GE. PSEL  .AND.  PDSE1L .LE. PSER) THEN
            DSE1L=DSETAN(1)
         ELSE
            DSE1L=FBCDSE(PDSE1L)
         END IF
         IF (PDSE1R .GE. PSEL  .AND.  PDSE1R .LE. PSER) THEN
            DSE1R=DSETAN(1)
         ELSE
            DSE1R=FBCDSE(PDSE1R)
         END IF
         DDSE1T(1)=(DSE1R - DSE1L)/(PDSE1R - PDSE1L)
         IF (PDSE2L .GE. PSEL  .AND.  PDSE2L .LE. PSER) THEN
            DSE2L=DSETAN(1)
         ELSE
            DSE2L=FBCDSE(PDSE2L)
         END IF
         IF (PDSE2R .GE. PSEL  .AND.  PDSE2R .LE. PSER) THEN
            DSE2R=DSETAN(1)
         ELSE
            DSE2R=FBCDSE(PDSE2R)
         END IF
         DDSE2T(1)=(DSE2R - DSE2L)/(PDSE2R - PDSE2L)
         KRL=FBCKR(PKRL)
         KRR=FBCKR(PKRR)
         DKRTAN=(KRR - KRL)/(PKRR - PKRL)
      END IF
C
      RETURN
      END
