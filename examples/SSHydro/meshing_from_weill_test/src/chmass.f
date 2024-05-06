C
C**************************  CHMASS ************************************
C
C  calculate soil moisture characteristics needed for 
C  mass balance calculations (general storage term)
C
C***********************************************************************
C
      SUBROUTINE CHMASS(NLKP,N,NT,NTRI,IVGHU,TETRA,TP,
     1                  PNEW,SNODI,PNODI,ETAI,ETAE)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I
      INTEGER  NLKP,N,NT,NTRI,IVGHU
      INTEGER  TETRA(5,*),TP(*)
      REAL*8   FVGSE,FVGDSE,FXVDMC,FHUSE,FHUDSE,FBCSE,FBCDSE
      REAL*8   PSI,SE,SWi
      REAL*8   PNEW(*),SNODI(*),PNODI(*),ETAI(*),ETAE(*)
      INCLUDE 'SOILCHAR.H'
C
      IF (IVGHU .EQ. -1) THEN
         CALL MOISTAB(NLKP,NT,NTRI,TETRA,PNEW,PNODI,SNODI,ETAE)
         CALL ELTNOD(N,NT,TP,TETRA,ETAE,ETAI)
      ELSE IF (IVGHU .EQ. 0) THEN
         DO I=1,N
            PSI      = PNEW(I)
            SE       = FVGSE(PSI,I)
            SWi      = VGPNOT(I)*SE + VGRMC(I)/PNODI(I)
            ETAI(I)  = SWi*SNODI(I) + 
     1                 PNODI(I)*VGPNOT(I)*FVGDSE(PSI,I)
         END DO 
      ELSE IF (IVGHU .EQ. 1) THEN
         DO I=1,N
            PSI      = PNEW(I)
            ETAI(I)  = FXVDMC(PSI,SNODI(I),PNODI(I),I)
         END DO
      ELSE IF (IVGHU .EQ. 2) THEN
         DO I=1,N
            PSI      = PNEW(I)
            SE       = FHUSE(PSI)
            SWi      = HUSWR1*SE + HUSWR
            ETAI(I)  = SWi*SNODI(I) + 
     1                 PNODI(I)*HUSWR1*FHUDSE(PSI)
         END DO
      ELSE IF (IVGHU .EQ. 3) THEN
         DO I=1,N
            PSI      = PNEW(I)
            SE       = FHUSE(PSI)
            SWi      = HUSWR1*SE + HUSWR
            ETAI(I)  = SWi*SNODI(I) + 
     1                 PNODI(I)*HUSWR1*FHUDSE(PSI)
         END DO
      ELSE 
         DO I=1,N
            PSI      = PNEW(I)
            SE       = FBCSE(PSI)
            SWi      = BCPORM(I)*SE + BCRMC/PNODI(I)
            ETAI(I)  = SWi*SNODI(I) + 
     1                 PNODI(I)*BCPORM(I)*FBCDSE(PSI)
         END DO 
      END IF
C
      RETURN
      END
