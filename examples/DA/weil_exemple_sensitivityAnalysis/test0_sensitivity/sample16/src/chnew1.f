C
C**************************  CHNEW1 ************************************
C
C  calculate soil moisture characteristics needed for Newton scheme
C  KSLOPE=1 : chord slope and analytical differentiation of moisture 
C             curves
C
C***********************************************************************
C
      SUBROUTINE CHNEW1(N,PTNEW,PTOLD,SNODI,PNODI,SW,CKRW,ETAI,
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
      REAL*8   DSENEW,DSEOLD,DMCNEW,DMCOLD
      REAL*8   PTNEW(*),PTOLD(*),SNODI(*),PNODI(*)
      REAL*8   SW(*),CKRW(*),ETAI(*),DCKRW(*),DETAI(*)
      INCLUDE 'SOILCHAR.H'
C
      IF (IVGHU .EQ. 0) THEN
         DO I=1,N
            PSI        = PTNEW(I)
            SE         = FVGSE(PSI,I)
            SW(I)      = VGPNOT(I)*SE + VGRMC(I)/PNODI(I)
            CKRW(I)    = FVGKR(PSI,SE,I)
            POLD       = PTOLD(I)
            DP         = PSI - POLD
            IF (DABS(DP) .LT. TOLKSL) THEN
               DSEDP   = FVGDSE(PSI,I)
               DSWDP   = VGPNOT(I)*DSEDP
               ETAI(I) = SW(I)*SNODI(I) + PNODI(I)*DSWDP
               DETAI(I)= DSWDP*SNODI(I) + 
     1                   PNODI(I)*VGPNOT(I)*FVGDDS(PSI,I)
               DCKRW(I)= FVGDKR(PSI,I)
            ELSE
               SEOLD   = FVGSE(POLD,I)
               DSENEW  = FVGDSE(PSI,I)
               DSEOLD  = FVGDSE(POLD,I)
               DSWDP   = VGPNOT(I)*((SE - SEOLD)/DP)
               KROLD   = FVGKR(POLD,SEOLD,I)
               ETAI(I) = SW(I)*SNODI(I) + PNODI(I)*DSWDP
               DETAI(I)= DSWDP*SNODI(I) + 
     1                   PNODI(I)*VGPNOT(I)*((DSENEW-DSEOLD)/DP)
               DCKRW(I)= (CKRW(I) - KROLD)/DP
            END IF
         END DO
      ELSE IF (IVGHU .EQ. 1) THEN
         DO I=1,N
            PSI        = PTNEW(I)
            SW(I)      = FXVMC(PSI,SNODI(I),PNODI(I),I)/PNODI(I)
            CKRW(I)    = FXVKR(PSI,I)
            POLD       = PTOLD(I)
            DP         = PSI - POLD
            IF (DABS(DP) .LT. TOLKSL) THEN
               ETAI(I) = FXVDMC(PSI,SNODI(I),PNODI(I),I)
               DETAI(I)= FXVDDM(PSI,PNODI(I),I)
               DCKRW(I)= FXVDKR(PSI,I)
            ELSE
               MC      = FXVMC(PSI,SNODI(I),PNODI(I),I)
               MCOLD   = FXVMC(POLD,SNODI(I),PNODI(I),I)
               DMCNEW  = FXVDMC(PSI,SNODI(I),PNODI(I),I)
               DMCOLD  = FXVDMC(POLD,SNODI(I),PNODI(I),I)
               KROLD   = FXVKR(POLD,I)
               ETAI(I) = (MC-MCOLD)/DP
               DETAI(I)= (DMCNEW-DMCOLD)/DP
               DCKRW(I)= (CKRW(I) - KROLD)/DP
            END IF
         END DO
      ELSE IF (IVGHU .EQ. 2) THEN
         DO I=1,N
            PSI        = PTNEW(I)
            SE         = FHUSE(PSI)
            SW(I)      = HUSWR1*SE + HUSWR
            CKRW(I)    = FHUKR2(PSI,SE)
            POLD       = PTOLD(I)
            DP         = PSI - POLD
            IF (DABS(DP) .LT. TOLKSL) THEN
               DSEDP   = FHUDSE(PSI)
               DSWDP   = HUSWR1*DSEDP
               ETAI(I) = SW(I)*SNODI(I) + PNODI(I)*DSWDP
               DETAI(I)= DSWDP*SNODI(I) + 
     1                   PNODI(I)*HUSWR1*FHUDDS(PSI)
               DCKRW(I)= FHUDK2(PSI,SE,DSEDP)
            ELSE
               SEOLD   = FHUSE(POLD)
               DSENEW  = FHUDSE(PSI)
               DSEOLD  = FHUDSE(POLD)
               DSWDP   = HUSWR1*((SE - SEOLD)/DP)
               KROLD   = FHUKR2(POLD,SEOLD)
               ETAI(I) = SW(I)*SNODI(I) + PNODI(I)*DSWDP
               DETAI(I)= DSWDP*SNODI(I) + 
     1                   PNODI(I)*HUSWR1*((DSENEW-DSEOLD)/DP)
               DCKRW(I)= (CKRW(I) - KROLD)/DP
            END IF
         END DO
      ELSE IF (IVGHU .EQ. 3) THEN
         DO I=1,N
            PSI        = PTNEW(I)
            SE         = FHUSE(PSI)
            SW(I)      = HUSWR1*SE + HUSWR
            CKRW(I)    = FHUKR3(PSI,SE)
            POLD       = PTOLD(I)
            DP         = PSI - POLD
            IF (DABS(DP) .LT. TOLKSL) THEN
               DSEDP   = FHUDSE(PSI)
               DSWDP   = HUSWR1*DSEDP
               ETAI(I) = SW(I)*SNODI(I) + PNODI(I)*DSWDP
               DETAI(I)= DSWDP*SNODI(I) + 
     1                   PNODI(I)*HUSWR1*FHUDDS(PSI)
               DCKRW(I)= FHUDK3(PSI,SE,DSEDP,CKRW(I))
            ELSE
               SEOLD   = FHUSE(POLD)
               DSENEW  = FHUDSE(PSI)
               DSEOLD  = FHUDSE(POLD)
               DSWDP   = HUSWR1*((SE - SEOLD)/DP)
               KROLD   = FHUKR3(POLD,SEOLD)
               ETAI(I) = SW(I)*SNODI(I) + PNODI(I)*DSWDP
               DETAI(I)= DSWDP*SNODI(I) + 
     1                   PNODI(I)*HUSWR1*((DSENEW-DSEOLD)/DP)
               DCKRW(I)= (CKRW(I) - KROLD)/DP
            END IF
         END DO
      ELSE
         DO I=1,N
            PSI        = PTNEW(I)
            SE         = FBCSE(PSI)
            SW(I)      = BCPORM(I)*SE + BCRMC/PNODI(I)
            CKRW(I)    = FBCKR(PSI)
            POLD       = PTOLD(I)
            DP         = PSI - POLD
            IF (DABS(DP) .LT. TOLKSL) THEN
               DSEDP   = FBCDSE(PSI)
               DSWDP   = BCPORM(I)*DSEDP
               ETAI(I) = SW(I)*SNODI(I) + PNODI(I)*DSWDP
               DETAI(I)= DSWDP*SNODI(I) + 
     1                   PNODI(I)*BCPORM(I)*FBCDDS(PSI)
               DCKRW(I)= FBCDKR(PSI)
            ELSE
               SEOLD   = FBCSE(POLD)
               DSENEW  = FBCDSE(PSI)
               DSEOLD  = FBCDSE(POLD)
               DSWDP   = BCPORM(I)*((SE - SEOLD)/DP)
               KROLD   = FBCKR(POLD)
               ETAI(I) = SW(I)*SNODI(I) + PNODI(I)*DSWDP
               DETAI(I)= DSWDP*SNODI(I) + 
     1                   PNODI(I)*BCPORM(I)*((DSENEW-DSEOLD)/DP)
               DCKRW(I)= (CKRW(I) - KROLD)/DP
            END IF
         END DO
      END IF
C
      RETURN
      END
