C
C**************************  PEATCH ************************************
C
C  Calculate soil moisture characteristics and porosity needed for the
C  Picard scheme for the case of a deformable peat soil.
C  Only the KSLOPE = 0 case (analytical differentiation of moisture
C  curves) is currently implemented.
C  Only IVGHU = -1, 0, or 4 soil moisture curve options are currently
C  implemented.
C
C***********************************************************************
C
      SUBROUTINE PEATCH(NLKP,N,NNOD,NT,NTRI,IVGHU,TETRA,TP,
     1                  PTNEW,Z,DEF,INDE,INDE0,PORE,
     2                  SNODI,PNODI,SW,CKRW,SENODI,ETAI,
     3                  CKRWE,SEELT,ETAE)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I,NS
      INTEGER  NLKP,N,NNOD,NT,NTRI,IVGHU
      INTEGER  TETRA(5,*),TP(*)
      REAL*8   FCINDE,FCDPORE
      REAL*8   FVGKR,FVGSE,FVGDSE,FBCKR,FBCSE,FBCDSE
      REAL*8   PSI,SE,D,CBETA
      REAL*8   PTNEW(*),Z(*),DEF(*),INDE(*),INDE0(*),PORE(*)
      REAL*8   SNODI(*),PNODI(*),SW(*),CKRW(*),SENODI(*),ETAI(*)
      REAL*8   CKRWE(*),SEELT(*),ETAE(*)
      INCLUDE 'SOILCHAR.H'
C
      IF (IVGHU .EQ. -1) THEN
         CALL MOISTABPICPT(NLKP,NT,NTRI,TETRA,PTNEW,SEELT,CKRWE,ETAE)
         CALL ELTNOD(N,NT,TP,TETRA,ETAE,ETAI)
         CALL ELTNOD(N,NT,TP,TETRA,CKRWE,CKRW)
         CALL ELTNOD(N,NT,TP,TETRA,SEELT,SENODI)
         DO I=1,N
            NS       = I-(NNOD*INT((I-1)/NNOD))
            D        = Z(NS) - Z(I)
            CBETA    = CBETA0+CANG*D
            IF (CBETA .GT. 1.0D0) CBETA=1.0D0
            PSI      = PTNEW(I)
            INDE(I)  = FCINDE(INDE(I),SW(I),D,PSI,PNODI(I),SNODI(I),
     1                 PORE(I))
            PORE(I)  = INDE(I)/(1.0D0 + INDE(I))
            VGPNOT(I)= (PORE(I)-VGRMC(I))/PORE(I)
            SE       = SENODI(I)
            SW(I)    = VGPNOT(I)*SE + VGRMC(I)/PORE(I)
            ETAI(I)  = SW(I)*SNODI(I) + (SW(I)*
     1                 FCDPORE(INDE(I),SW(I),D,PNODI(I)) + 1.0D0)*
     2                 ETAI(I)*PORE(I)*VGPNOT(I)
            IF (SW(I) .LT. 1.0D0) ETAI(I) = ETAI(I) - SW(I)*SNODI(I)
            DEF(I)   = ((1.0D0+INDE(I))/(1.0D0+INDE0(I)))**CBETA
         END DO
         CALL NODELT(NT,TETRA,ETAI,ETAE)
      ELSE IF (IVGHU .EQ. 0) THEN
         DO I=1,N
            NS       = I-(NNOD*INT((I-1)/NNOD))
            D        = Z(NS) - Z(I)
            CBETA    = CBETA0+CANG*D
            IF (CBETA .GT. 1.0D0) CBETA=1.0D0
            PSI      = PTNEW(I)
            SE       = FVGSE(PSI,I)
            INDE(I)  = FCINDE(INDE(I),SW(I),D,PSI,PNODI(I),SNODI(I),
     1                 PORE(I))
            PORE(I)  = INDE(I)/(1.0D0 + INDE(I))
            VGPNOT(I)= (PORE(I)-VGRMC(I))/PORE(I)
            SW(I)    = VGPNOT(I)*SE + VGRMC(I)/PORE(I)
            ETAI(I)  = SW(I)*SNODI(I) + (SW(I)*
     1                 FCDPORE(INDE(I),SW(I),D,PNODI(I)) + 1.0D0)*
     2                 FVGDSE(PSI,I)*PORE(I)*VGPNOT(I)
            IF (SW(I) .LT. 1.0D0) ETAI(I) = ETAI(I) - SW(I)*SNODI(I)
            CKRW(I)  = FVGKR(PSI,SE,I)
            DEF(I)   = ((1.0D0+INDE(I))/(1.0D0+INDE0(I)))**CBETA
         END DO
         CALL NODELT(NT,TETRA,CKRW,CKRWE)
         CALL NODELT(NT,TETRA,ETAI,ETAE)
      ELSE IF(IVGHU .EQ. 4) THEN
         DO I=1,N
            NS       = I-(NNOD*INT((I-1)/NNOD))
            D        = Z(NS) - Z(I)
            CBETA    = CBETA0+CANG*D
            IF (CBETA .GT. 1.0D0) CBETA=1.0D0
            PSI      = PTNEW(I)
            SE       = FBCSE(PSI)
            INDE(I)  = FCINDE(INDE(I),SW(I),D,PSI,PNODI(I),SNODI(I),
     1                 PORE(I)) 
            PORE(I)  = INDE(I)/(1.0D0 + INDE(I))
            SW(I)    = BCPORM(I)*SE + BCRMC/PORE(I)
            ETAI(I)  = SW(I)*SNODI(I) + (SW(I)*
     1                 FCDPORE(INDE(I),SW(I),D,PNODI(I)) + 1.0D0)*
     2                 FBCDSE(PSI)*PORE(I)*BCPORM(I)
            CKRW(I)  = FBCKR(PSI)
            DEF(I)   = ((1.0D0+INDE(I))/(1.0D0+INDE0(I)))**CBETA
         END DO
         CALL NODELT(NT,TETRA,CKRW,CKRWE)
         CALL NODELT(NT,TETRA,ETAI,ETAE)
      END IF
C
      RETURN
      END
