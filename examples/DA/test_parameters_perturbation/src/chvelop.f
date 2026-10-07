C
C**************************  CHVELOP ***********************************
C
C  calculate soil moisture characteristics needed for storage and
C  velocity calculations and for output (peat soil case)
C
C***********************************************************************
C
      SUBROUTINE CHVELOP(NLKP,NNOD,N,NT,NTRI,IVGHU,TETRA,TP,
     1                   PNEW,Z,PORE,SW,CKRW,SENODI,CKRWE,SEELT)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I,NS
      INTEGER  NLKP,NNOD,N,NT,NTRI,IVGHU
      INTEGER  TETRA(5,*),TP(*)
      REAL*8   FVGKR,FVGSE,FBCKR,FBCSE
      REAL*8   PSI,SE,D,CBETA
      REAL*8   PNEW(*),Z(*),PORE(*)
      REAL*8   SW(*),CKRW(*),SENODI(*),CKRWE(*),SEELT(*)
      INCLUDE 'SOILCHAR.H'
C
      IF (IVGHU .EQ. -1) THEN
         CALL MOISTABVELPT(NLKP,NT,NTRI,TETRA,PNEW,SEELT,CKRWE)
         CALL ELTNOD(N,NT,TP,TETRA,CKRWE,CKRW)
         CALL ELTNOD(N,NT,TP,TETRA,SEELT,SENODI)
         DO I=1,N
            NS       = I-(NNOD*INT((I-1)/NNOD))
            D        = Z(NS) - Z(I)
            CBETA    = CBETA0+CANG*D
            IF (CBETA .GT. 1.0D0) CBETA=1.0D0
            PSI      = PNEW(I)
            SE       = SENODI(I)
            VGPNOT(I)= (PORE(I)-VGRMC(I))/PORE(I)
            SW(I)    = VGPNOT(I)*SE + VGRMC(I)/PORE(I)
         END DO 
      ELSE IF (IVGHU .EQ. 0) THEN
         DO I=1,N
            NS       = I-(NNOD*INT((I-1)/NNOD))
            D        = Z(NS) - Z(I)
            CBETA    = CBETA0+CANG*D
            IF (CBETA .GT. 1.0D0) CBETA=1.0D0
            PSI      = PNEW(I)
            SE       = FVGSE(PSI,I)
            VGPNOT(I)= (PORE(I)-VGRMC(I))/PORE(I)
            SW(I)    = VGPNOT(I)*SE + VGRMC(I)/PORE(I)
            CKRW(I)  = FVGKR(PSI,SE,I)
          END DO
      ELSE IF (IVGHU .EQ. 4) THEN
         DO I=1,N
            NS       = I-(NNOD*INT((I-1)/NNOD))
            D        = Z(NS) - Z(I)
            CBETA    = CBETA0+CANG*D
            IF (CBETA .GT. 1.0D0) CBETA=1.0D0
            PSI      = PNEW(I)
            SE       = FBCSE(PSI)
            SW(I)    = BCPORM(I)*SE + BCRMC/PORE(I)
            CKRW(I)  = FBCKR(PSI)
         END DO 
      END IF
C
      RETURN
      END
