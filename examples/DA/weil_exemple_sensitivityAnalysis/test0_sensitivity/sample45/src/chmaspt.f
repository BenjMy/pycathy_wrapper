C
C**************************  CHMASPT ***********************************
C
C  calculate soil moisture characteristics needed for 
C  mass balance calculations (general storage term)
C  in case of deformable peat soil
C
C***********************************************************************
C
      SUBROUTINE CHMASPT(NLKP,N,NNOD,NT,NTRI,IVGHU,TETRA,TP,
     1                   PNEW,SNODI,PNODI,PORE,INDE,Z,
     2                   SENODI,ETAI,SEELT,ETAE)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I,NS
      INTEGER  NLKP,N,NNOD,NT,NTRI,IVGHU
      INTEGER  TETRA(5,*),TP(*)
      REAL*8   FCDPORE,FVGSE,FVGDSE,FBCSE,FBCDSE
      REAL*8   PSI,SE,SWi,D
      REAL*8   PNEW(*),SNODI(*),PNODI(*),PORE(*),INDE(*),Z(*)
      REAL*8   SENODI(*),ETAI(*),SEELT(*),ETAE(*)
      INCLUDE 'SOILCHAR.H'
C
      IF (IVGHU .EQ. -1) THEN
         CALL MOISTABPT(NLKP,NT,NTRI,TETRA,PNEW,SEELT,ETAE)             
         CALL ELTNOD(N,NT,TP,TETRA,ETAE,ETAI)
         CALL ELTNOD(N,NT,TP,TETRA,SEELT,SENODI)
         DO I=1,N
            NS       = I-(NNOD*INT((I-1)/NNOD))
            D        = Z(NS) - Z(I)
            PSI      = PNEW(I)
            SE       = SENODI(I)
            VGPNOT(I)= (PORE(I)-VGRMC(I))/PORE(I)
            SWi      = VGPNOT(I)*SE + VGRMC(I)/PORE(I)
            ETAI(I)  = SWi*SNODI(I) +
     1                 (SWi*FCDPORE(INDE(I),SWi,D,PNODI(I)) + 
     2                  PORE(I)*VGPNOT(I)) * ETAI(I)
            IF (SWi .EQ. 1.0D0) ETAI(I) = ETAI(I) - SWi*SNODI(I)
         END DO
      ELSE IF (IVGHU .EQ. 0) THEN
         DO I=1,N
            NS       = I-(NNOD*INT((I-1)/NNOD))
            D        = Z(NS) - Z(I)
            PSI      = PNEW(I)
            SE       = FVGSE(PSI,I)
            VGPNOT(I)= (PORE(I)-VGRMC(I))/PORE(I)
            SWi      = VGPNOT(I)*SE + VGRMC(I)/PORE(I) 
            ETAI(I)  = SWi*SNODI(I) +
     1                 (SWi*FCDPORE(INDE(I),SWi,D,PNODI(I)) +
     2                  PORE(I)*VGPNOT(I)) * FVGDSE(PSI,I) 
            IF (SWi .EQ. 1.0D0) ETAI(I) = ETAI(I) - SWi*SNODI(I)
         END DO 
      ELSE IF (IVGHU .EQ. 4) THEN
         DO I=1,N
            NS       = I-(NNOD*INT((I-1)/NNOD))
            D        = Z(NS) - Z(I)
            PSI      = PNEW(I)
            SE       = FBCSE(PSI)
            SWi      = BCPORM(I)*SE + BCRMC/PORE(I)
            ETAI(I)  = SWi*SNODI(I) +
     1                 (SWi*FCDPORE(INDE(I),SWi,D,PNODI(I)) +
     2                  PORE(I)*BCPORM(I)) * FBCDSE(PSI) 
         END DO 
      END IF
C
      RETURN
      END
