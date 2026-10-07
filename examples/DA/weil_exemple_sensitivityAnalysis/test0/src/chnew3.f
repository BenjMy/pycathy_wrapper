C
C**************************  CHNEW3 ************************************
C
C  calculate soil moisture characteristics needed for Newton scheme
C  KSLOPE=3 : localized chord slope and analytical differentiation of
C             moisture curves
C
C***********************************************************************
C
      SUBROUTINE CHNEW3(N,PTNEW,PTOLD,SNODI,PNODI,SW,CKRW,ETAI,
     1                  DCKRW,DETAI,IVGHU)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I
      INTEGER  N,IVGHU
      REAL*8   FVGKR,FVGDKR,FVGSE,FVGDSE,FVGDDS
      REAL*8   FXVKR,FXVDKR,FXVMC,FXVDMC,FXVDDM
      REAL*8   FHUKR2,FHUKR3,FHUDK2,FHUDK3,FHUSE,FHUDSE,FHUDDS
      REAL*8   FBCKR,FBCDKR,FBCSE,FBCDSE,FBCDDS
      REAL*8   PSI,POLD,DP,SE,SEOLD,MC,MCOLD,KROLD,DSEDP,DSWDP
      REAL*8   DMCNEW,DMCOLD
      REAL*8   PTNEW(*),PTOLD(*),SNODI(*),PNODI(*)
      REAL*8   SW(*),CKRW(*),ETAI(*),DCKRW(*),DETAI(*)
      INCLUDE 'SOILCHAR.H'
C
      IF (IVGHU .EQ. 0) THEN
         DO I=1,N
            PSI           = PTNEW(I)
            SE            = FVGSE(PSI,I)
            SW(I)         = VGPNOT(I)*SE + VGRMC(I)/PNODI(I)
            CKRW(I)       = FVGKR(PSI,SE,I)
            IF (PSI .GE. PSEL  .AND.  PSI .LE. PSER) THEN
               POLD       = PTOLD(I)
               DP         = PSI - POLD
               IF (DABS(DP) .LT. TOLKSL) THEN
                  DSWDP   = VGPNOT(I)*FVGDSE(PSI,I)
                  ETAI(I) = SW(I)*SNODI(I) + PNODI(I)*DSWDP
               ELSE
                  SEOLD   = FVGSE(POLD,I)
                  DSWDP   = VGPNOT(I)*((SE - SEOLD)/DP)
                  ETAI(I) = SW(I)*SNODI(I) + PNODI(I)*DSWDP
               END IF
            ELSE
               DSWDP      = VGPNOT(I)*FVGDSE(PSI,I)
               ETAI(I)    = SW(I)*SNODI(I) + PNODI(I)*DSWDP
            END IF
            IF ((PSI .GE. PDSE1L  .AND.  PSI .LE. PDSE1R) .OR.
     1          (PSI .GE. PDSE2L  .AND.  PSI .LE. PDSE2R)) THEN
               POLD       = PTOLD(I)
               DP         = PSI - POLD
               IF (DABS(DP) .LT. TOLKSL) THEN
                  DETAI(I)= DSWDP*SNODI(I) + 
     1                      PNODI(I)*VGPNOT(I)*FVGDDS(PSI,I)
               ELSE
                  DETAI(I)= DSWDP*SNODI(I) + PNODI(I)*
     1                     VGPNOT(I)*((FVGDSE(PSI,I)-FVGDSE(POLD,I))/DP)
               END IF
            ELSE
               DETAI(I)   = DSWDP*SNODI(I) + 
     1                      PNODI(I)*VGPNOT(I)*FVGDDS(PSI,I)
            END IF
            IF (PSI .GE. PKRL  .AND.  PSI .LE. PKRR) THEN
               POLD       = PTOLD(I)
               DP         = PSI - POLD
               IF (DABS(DP) .LT. TOLKSL) THEN
                  DCKRW(I)= FVGDKR(PSI,I)
               ELSE
                  SEOLD   = FVGSE(POLD,I)
                  DCKRW(I)= (CKRW(I) - FVGKR(POLD,SEOLD,I))/DP
               END IF
            ELSE
               DCKRW(I)   = FVGDKR(PSI,I)
            END IF
         END DO
      ELSE IF (IVGHU .EQ. 1) THEN
         DO I=1,N
            PSI           = PTNEW(I)
            SW(I)         = FXVMC(PSI,SNODI(I),PNODI(I),I)/PNODI(I)
            CKRW(I)       = FXVKR(PSI,I)
            IF (PSI .GE. PSEL  .AND.  PSI .LE. PSER) THEN
               POLD       = PTOLD(I)
               DP         = PSI - POLD
               IF (DABS(DP) .LT. TOLKSL) THEN
                  ETAI(I) = FXVDMC(PSI,SNODI(I),PNODI(I),I)
               ELSE
                  MC      = FXVMC(PSI,SNODI(I),PNODI(I),I)
                  MCOLD   = FXVMC(POLD,SNODI(I),PNODI(I),I)
                  ETAI(I) = (MC-MCOLD)/DP
               END IF 
            ELSE
               ETAI(I)    = FXVDMC(PSI,SNODI(I),PNODI(I),I)
            END IF
            IF ((PSI .GE. PDSE1L  .AND.  PSI .LE. PDSE1R) .OR.
     1          (PSI .GE. PDSE2L  .AND.  PSI .LE. PDSE2R)) THEN
               POLD       = PTOLD(I)
               DP         = PSI - POLD
               IF (DABS(DP) .LT. TOLKSL) THEN
                  DETAI(I)= FXVDDM(PSI,PNODI(I),I)
               ELSE
                  DMCNEW  = FXVDMC(PSI,SNODI(I),PNODI(I),I)
                  DMCOLD  = FXVDMC(POLD,SNODI(I),PNODI(I),I)
                  DETAI(I)= (DMCNEW-DMCOLD)/DP
               END IF
            ELSE
               DETAI(I)   = FXVDDM(PSI,PNODI(I),I)
            END IF
            IF (PSI .GE. PKRL  .AND.  PSI .LE. PKRR) THEN
               POLD       = PTOLD(I)
               DP         = PSI - POLD
               IF (DABS(DP) .LT. TOLKSL) THEN
                  DCKRW(I)= FXVDKR(PSI,I)
               ELSE
                  DCKRW(I)= (CKRW(I) - FXVKR(POLD,I))/DP
               END IF
            ELSE
               DCKRW(I)   = FXVDKR(PSI,I)
            END IF
         END DO
      ELSE IF (IVGHU .EQ. 2) THEN
         DO I=1,N
            PSI           = PTNEW(I)
            SE            = FHUSE(PSI)
            SW(I)         = HUSWR1*SE + HUSWR
            CKRW(I)       = FHUKR2(PSI,SE)
            IF (PSI .GE. PSEL  .AND.  PSI .LE. PSER) THEN
               POLD       = PTOLD(I)
               DP         = PSI - POLD
               IF (DABS(DP) .LT. TOLKSL) THEN
                  DSEDP   = FHUDSE(PSI)
                  DSWDP   = HUSWR1*DSEDP
                  ETAI(I) = SW(I)*SNODI(I) + PNODI(I)*DSWDP
               ELSE
                  SEOLD   = FHUSE(POLD)
                  DSEDP   = (SE - SEOLD)/DP
                  DSWDP   = HUSWR1*DSEDP
                  ETAI(I) = SW(I)*SNODI(I) + PNODI(I)*DSWDP
               END IF
            ELSE
               DSEDP      = FHUDSE(PSI)
               DSWDP      = HUSWR1*DSEDP
               ETAI(I)    = SW(I)*SNODI(I) + PNODI(I)*DSWDP
            END IF
            IF ((PSI .GE. PDSE1L  .AND.  PSI .LE. PDSE1R) .OR.
     1          (PSI .GE. PDSE2L  .AND.  PSI .LE. PDSE2R)) THEN
               POLD       = PTOLD(I)
               DP         = PSI - POLD
               IF (DABS(DP) .LT. TOLKSL) THEN
                  DETAI(I)= DSWDP*SNODI(I) + 
     1                      PNODI(I)*HUSWR1*FHUDDS(PSI)
               ELSE
                  DETAI(I)= DSWDP*SNODI(I) + PNODI(I)*
     1                      HUSWR1*((FHUDSE(PSI)-FHUDSE(POLD))/DP)
               END IF
            ELSE
               DETAI(I)   = DSWDP*SNODI(I) + 
     1                      PNODI(I)*HUSWR1*FHUDDS(PSI)
            END IF
            IF (PSI .GE. PKRL  .AND.  PSI .LE. PKRR) THEN
               POLD       = PTOLD(I)
               DP         = PSI - POLD
               IF (DABS(DP) .LT. TOLKSL) THEN
                  DCKRW(I)= FHUDK2(PSI,SE,DSEDP)
               ELSE
                  KROLD   = FHUKR2(POLD,SEOLD)
                  DCKRW(I)= (CKRW(I) - KROLD)/DP
               END IF
            ELSE
               DCKRW(I)   = FHUDK2(PSI,SE,DSEDP)
            END IF
         END DO
      ELSE IF (IVGHU .EQ. 3) THEN
         DO I=1,N
            PSI           = PTNEW(I)
            SE            = FHUSE(PSI)
            SW(I)         = HUSWR1*SE + HUSWR
            CKRW(I)       = FHUKR3(PSI,SE)
            IF (PSI .GE. PSEL  .AND.  PSI .LE. PSER) THEN
               POLD       = PTOLD(I)
               DP         = PSI - POLD
               IF (DABS(DP) .LT. TOLKSL) THEN
                  DSEDP   = FHUDSE(PSI)
                  DSWDP   = HUSWR1*DSEDP
                  ETAI(I) = SW(I)*SNODI(I) + PNODI(I)*DSWDP
               ELSE
                  SEOLD   = FHUSE(POLD)
                  DSEDP   = (SE - SEOLD)/DP
                  DSWDP   = HUSWR1*DSEDP
                  ETAI(I) = SW(I)*SNODI(I) + PNODI(I)*DSWDP
               END IF
            ELSE
               DSEDP      = FHUDSE(PSI)
               DSWDP      = HUSWR1*DSEDP
               ETAI(I)    = SW(I)*SNODI(I) + PNODI(I)*DSWDP
            END IF
            IF ((PSI .GE. PDSE1L  .AND.  PSI .LE. PDSE1R) .OR.
     1          (PSI .GE. PDSE2L  .AND.  PSI .LE. PDSE2R)) THEN
               POLD       = PTOLD(I)
               DP         = PSI - POLD
               IF (DABS(DP) .LT. TOLKSL) THEN
                  DETAI(I)= DSWDP*SNODI(I) + 
     1                      PNODI(I)*HUSWR1*FHUDDS(PSI)
               ELSE
                  DETAI(I)= DSWDP*SNODI(I) + PNODI(I)*
     1                      HUSWR1*((FHUDSE(PSI)-FHUDSE(POLD))/DP)
               END IF
            ELSE
               DETAI(I)   = DSWDP*SNODI(I) + 
     1                      PNODI(I)*HUSWR1*FHUDDS(PSI)
            END IF
            IF (PSI .GE. PKRL  .AND.  PSI .LE. PKRR) THEN
               POLD       = PTOLD(I)
               DP         = PSI - POLD
               IF (DABS(DP) .LT. TOLKSL) THEN
                  DCKRW(I)= FHUDK3(PSI,SE,DSEDP,CKRW(I))
               ELSE
                  KROLD   = FHUKR3(POLD,SEOLD)
                  DCKRW(I)= (CKRW(I) - KROLD)/DP
               END IF
            ELSE 
               DCKRW(I)   = FHUDK3(PSI,SE,DSEDP,CKRW(I))
            END IF
         END DO
      ELSE
         DO I=1,N
            PSI           = PTNEW(I)
            SE            = FBCSE(PSI)
            SW(I)         = BCPORM(I)*SE + BCRMC/PNODI(I)
            CKRW(I)       = FBCKR(PSI)
            IF (PSI .GE. PSEL  .AND.  PSI .LE. PSER) THEN
               POLD       = PTOLD(I)
               DP         = PSI - POLD
               IF (DABS(DP) .LT. TOLKSL) THEN
                  DSWDP   = BCPORM(I)*FBCDSE(PSI)
                  ETAI(I) = SW(I)*SNODI(I) + PNODI(I)*DSWDP
               ELSE
                  SEOLD   = FBCSE(POLD)
                  DSWDP   = BCPORM(I)*((SE - SEOLD)/DP)
                  ETAI(I) = SW(I)*SNODI(I) + PNODI(I)*DSWDP
               END IF
            ELSE
               DSWDP      = BCPORM(I)*FBCDSE(PSI)
               ETAI(I)    = SW(I)*SNODI(I) + PNODI(I)*DSWDP
            END IF
            IF ((PSI .GE. PDSE1L  .AND.  PSI .LE. PDSE1R) .OR.
     1          (PSI .GE. PDSE2L  .AND.  PSI .LE. PDSE2R)) THEN
               POLD       = PTOLD(I)
               DP         = PSI - POLD
               IF (DABS(DP) .LT. TOLKSL) THEN
                  DETAI(I)= DSWDP*SNODI(I) + 
     1                      PNODI(I)*BCPORM(I)*FBCDDS(PSI)
               ELSE
                  DETAI(I)= DSWDP*SNODI(I) + PNODI(I)*
     1                      BCPORM(I)*((FBCDSE(PSI)-FBCDSE(POLD))/DP)
               END IF
            ELSE
               DETAI(I)   = DSWDP*SNODI(I) + 
     1                      PNODI(I)*BCPORM(I)*FBCDDS(PSI)
            END IF
            IF (PSI .GE. PKRL  .AND.  PSI .LE. PKRR) THEN
               POLD       = PTOLD(I)
               DP         = PSI - POLD
               IF (DABS(DP) .LT. TOLKSL) THEN
                  DCKRW(I)= FBCDKR(PSI)
               ELSE
                  SEOLD   = FBCSE(POLD)
                  DCKRW(I)= (CKRW(I) - FBCKR(POLD))/DP
               END IF
            ELSE
               DCKRW(I)   = FBCDKR(PSI)
            END IF
         END DO
      END IF
C
      RETURN
      END
