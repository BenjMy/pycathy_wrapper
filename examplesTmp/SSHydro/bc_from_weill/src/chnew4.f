C
C**************************  CHNEW4 ************************************
C
C  calculate soil moisture characteristics needed for Newton scheme
C  KSLOPE=4 : localized tangent slope differentiation of moisture 
C  curves. 
C  Note that DSETAN, DDSE1T, and DDSE2T contain tangent slope
C  values at each node only for the case IVGHU=1; for the other 
C  IVGHU cases the tangent slope values are constant for all nodes
C  and we use DSETAN(1), DDSE1T(1), and DDSE2T(1).
C
C***********************************************************************
C
      SUBROUTINE CHNEW4(N,PTNEW,SNODI,PNODI,SW,CKRW,ETAI,DCKRW,DETAI,
     1                  IVGHU)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I
      INTEGER  N,IVGHU
      REAL*8   FVGKR,FVGDKR,FVGSE,FVGDSE,FVGDDS
      REAL*8   FXVKR,FXVMC,FXVDKR,FXVDMC,FXVDDM
      REAL*8   FHUKR2,FHUKR3,FHUDK2,FHUDK3,FHUSE,FHUDSE,FHUDDS
      REAL*8   FBCKR,FBCDKR,FBCSE,FBCDSE,FBCDDS
      REAL*8   PSI,SE,DSEDP,DSWDP
      REAL*8   PTNEW(*),SNODI(*),PNODI(*)
      REAL*8   SW(*),CKRW(*),ETAI(*),DCKRW(*),DETAI(*)
      INCLUDE 'SOILCHAR.H'
C
      IF (IVGHU .EQ. 0) THEN
         DO I=1,N
            PSI        = PTNEW(I)
            SE         = FVGSE(PSI,I)
            SW(I)      = VGPNOT(I)*SE + VGRMC(I)/PNODI(I)
            CKRW(I)    = FVGKR(PSI,SE,I)
            IF (PSI .GE. PSEL  .AND.  PSI .LE. PSER) THEN
               DSWDP   = VGPNOT(I)*DSETAN(1)
               ETAI(I) = SW(I)*SNODI(I) + PNODI(I)*DSWDP
            ELSE
               DSWDP   = VGPNOT(I)*FVGDSE(PSI,I)
               ETAI(I) = SW(I)*SNODI(I) + PNODI(I)*DSWDP
            END IF
            IF (PSI .GE. PDSE1L  .AND.  PSI .LE. PDSE1R) THEN
               DETAI(I)= DSWDP*SNODI(I) + 
     1                   PNODI(I)*VGPNOT(I)*DDSE1T(1)
            ELSE IF (PSI .GE. PDSE2L  .AND.  PSI .LE. PDSE2R) THEN
               DETAI(I)= DSWDP*SNODI(I) + 
     1                   PNODI(I)*VGPNOT(I)*DDSE2T(1)
            ELSE
               DETAI(I)= DSWDP*SNODI(I) + 
     1                   PNODI(I)*VGPNOT(I)*FVGDDS(PSI,I)
            END IF
            IF (PSI .GE. PKRL  .AND.  PSI .LE. PKRR) THEN
               DCKRW(I)= DKRTAN
            ELSE
               DCKRW(I)= FVGDKR(PSI,I)
            END IF
         END DO
      ELSE IF (IVGHU .EQ. 1) THEN
         DO I=1,N
            PSI        = PTNEW(I)
            SW(I)      = FXVMC(PSI,SNODI(I),PNODI(I),I)/PNODI(I)
            CKRW(I)    = FXVKR(PSI,I)
            IF (PSI .GE. PSEL  .AND.  PSI .LE. PSER) THEN
               ETAI(I) = DSETAN(I)
            ELSE
               ETAI(I) = FXVDMC(PSI,SNODI(I),PNODI(I),I)
            END IF
            IF (PSI .GE. PDSE1L  .AND.  PSI .LE. PDSE1R) THEN
               DETAI(I)= DDSE1T(I)
            ELSE IF (PSI .GE. PDSE2L  .AND.  PSI .LE. PDSE2R) THEN
               DETAI(I)= DDSE2T(I)
            ELSE
               DETAI(I)= FXVDDM(PSI,PNODI(I),I)
            END IF
            IF (PSI .GE. PKRL  .AND.  PSI .LE. PKRR) THEN
               DCKRW(I)= DKRTAN
            ELSE
               DCKRW(I)= FXVDKR(PSI,I)
            END IF
         END DO
      ELSE IF (IVGHU .EQ. 2) THEN
         DO I=1,N
            PSI        = PTNEW(I)
            SE         = FHUSE(PSI)
            SW(I)      = HUSWR1*SE + HUSWR
            CKRW(I)    = FHUKR2(PSI,SE)
            IF (PSI .GE. PSEL  .AND.  PSI .LE. PSER) THEN
               DSEDP   = DSETAN(1)
               DSWDP   = HUSWR1*DSEDP
               ETAI(I) = SW(I)*SNODI(I) + PNODI(I)*DSWDP
            ELSE
               DSEDP   = FHUDSE(PSI)
               DSWDP   = HUSWR1*DSEDP
               ETAI(I) = SW(I)*SNODI(I) + PNODI(I)*DSWDP
            END IF
            IF (PSI .GE. PDSE1L  .AND.  PSI .LE. PDSE1R) THEN
               DETAI(I)= DSWDP*SNODI(I) + 
     1                   PNODI(I)*HUSWR1*DDSE1T(1)
            ELSE IF (PSI .GE. PDSE2L  .AND.  PSI .LE. PDSE2R) THEN
               DETAI(I)= DSWDP*SNODI(I) + 
     1                   PNODI(I)*HUSWR1*DDSE2T(1)
            ELSE
               DETAI(I)= DSWDP*SNODI(I) + 
     1                   PNODI(I)*HUSWR1*FHUDDS(PSI)
            END IF
            IF (PSI .GE. PKRL  .AND.  PSI .LE. PKRR) THEN
               DCKRW(I)= DKRTAN
            ELSE
               DCKRW(I)= FHUDK2(PSI,SE,DSEDP)
            END IF
         END DO
      ELSE IF (IVGHU .EQ. 3) THEN
         DO I=1,N
            PSI        = PTNEW(I)
            SE         = FHUSE(PSI)
            SW(I)      = HUSWR1*SE + HUSWR
            CKRW(I)    = FHUKR3(PSI,SE)
            IF (PSI .GE. PSEL  .AND.  PSI .LE. PSER) THEN
               DSEDP   = DSETAN(1)
               DSWDP   = HUSWR1*DSEDP
               ETAI(I) = SW(I)*SNODI(I) + PNODI(I)*DSWDP
            ELSE
               DSEDP   = FHUDSE(PSI)
               DSWDP   = HUSWR1*DSEDP
               ETAI(I) = SW(I)*SNODI(I) + PNODI(I)*DSWDP
            END IF
            IF (PSI .GE. PDSE1L  .AND.  PSI .LE. PDSE1R) THEN
               DETAI(I)= DSWDP*SNODI(I) + 
     1                   PNODI(I)*HUSWR1*DDSE1T(1)
            ELSE IF (PSI .GE. PDSE2L  .AND.  PSI .LE. PDSE2R) THEN
               DETAI(I)= DSWDP*SNODI(I) + 
     1                   PNODI(I)*HUSWR1*DDSE2T(1)
            ELSE
               DETAI(I)= DSWDP*SNODI(I) + 
     1                   PNODI(I)*HUSWR1*FHUDDS(PSI)
            END IF
            IF (PSI .GE. PKRL  .AND.  PSI .LE. PKRR) THEN
               DCKRW(I)= DKRTAN
            ELSE 
               DCKRW(I)= FHUDK3(PSI,SE,DSEDP,CKRW(I))
            END IF
         END DO
      ELSE
         DO I=1,N
            PSI        = PTNEW(I)
            SE         = FBCSE(PSI)
            SW(I)      = BCPORM(I)*SE + BCRMC/PNODI(I)
            CKRW(I)    = FBCKR(PSI)
            IF (PSI .GE. PSEL  .AND.  PSI .LE. PSER) THEN
               DSWDP   = BCPORM(I)*DSETAN(1)
               ETAI(I) = SW(I)*SNODI(I) + PNODI(I)*DSWDP
            ELSE
               DSWDP   = BCPORM(I)*FBCDSE(PSI)
               ETAI(I) = SW(I)*SNODI(I) + PNODI(I)*DSWDP
            END IF
            IF (PSI .GE. PDSE1L  .AND.  PSI .LE. PDSE1R) THEN
               DETAI(I)= DSWDP*SNODI(I) + 
     1                   PNODI(I)*BCPORM(I)*DDSE1T(1)
            ELSE IF (PSI .GE. PDSE2L  .AND.  PSI .LE. PDSE2R) THEN
               DETAI(I)= DSWDP*SNODI(I) + 
     1                   PNODI(I)*BCPORM(I)*DDSE2T(1)
            ELSE
               DETAI(I)= DSWDP*SNODI(I) + 
     1                   PNODI(I)*BCPORM(I)*FBCDDS(PSI)
            END IF
            IF (PSI .GE. PKRL  .AND.  PSI .LE. PKRR) THEN
               DCKRW(I)= DKRTAN
            ELSE
               DCKRW(I)= FBCDKR(PSI)
            END IF
         END DO
      END IF
C
      RETURN
      END
