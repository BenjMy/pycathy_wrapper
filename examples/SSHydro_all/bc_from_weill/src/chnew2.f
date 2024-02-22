C
C**************************  CHNEW2 ************************************
C
C  calculate soil moisture characteristics needed for Newton scheme
C  KSLOPE=2 : chord slope and centered difference formulas for 
C             differentiation of moisture curves
C
C***********************************************************************
C
      SUBROUTINE CHNEW2(N,PTNEW,PTOLD,SNODI,PNODI,SW,CKRW,ETAI,
     1                  DCKRW,DETAI,IVGHU)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I
      INTEGER  N,IVGHU
      REAL*8   FVGKR,FVGSE,FVGDSE,FXVKR,FXVMC,FXVDMC
      REAL*8   FHUKR2,FHUKR3,FHUSE,FHUDSE,FBCKR,FBCSE,FBCDSE
      REAL*8   PSI,POLD,PPDEL,PMDEL,DP,SE,SEOLD,SEPDEL,SEMDEL
      REAL*8   DSE,DDSE,DSWDP
      REAL*8   MC,MCOLD,MCPDEL,MCMDEL,DMCNEW,DMCOLD,TOL2
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
               TOL2    = 2.0D0*TOLKSL
               PPDEL   = PSI + TOLKSL
               PMDEL   = PSI - TOLKSL
               SEPDEL  = FVGSE(PPDEL,I)
               SEMDEL  = FVGSE(PMDEL,I)
               DSE     = (SEPDEL - SEMDEL)/TOL2
               DDSE    = (SEPDEL - 2.0D0*SE + SEMDEL)/(TOLKSL*TOLKSL)
               DCKRW(I)= (FVGKR(PPDEL,SEPDEL,I) - 
     1                   FVGKR(PMDEL,SEMDEL,I))/TOL2
            ELSE
               SEOLD   = FVGSE(POLD,I)
               DSE     = (SE - SEOLD)/DP
               DDSE    = (FVGDSE(PSI,I) - FVGDSE(POLD,I))/DP
               DCKRW(I)= (CKRW(I) - FVGKR(POLD,SEOLD,I))/DP
            END IF
            DSWDP      = VGPNOT(I)*DSE
            ETAI(I)    = SW(I)*SNODI(I) + PNODI(I)*DSWDP
            DETAI(I)   = DSWDP*SNODI(I) + PNODI(I)*VGPNOT(I)*DDSE
         END DO 
      ELSE IF (IVGHU .EQ. 1) THEN
         DO I=1,N
            PSI        = PTNEW(I)
            MC         = FXVMC(PSI,SNODI(I),PNODI(I),I)
            SW(I)      = MC/PNODI(I)
            CKRW(I)    = FXVKR(PSI,I)
            POLD       = PTOLD(I)
            DP         = PSI - POLD
            IF (DABS(DP) .LT. TOLKSL) THEN
               TOL2    = 2.0D0*TOLKSL
               PPDEL   = PSI + TOLKSL
               PMDEL   = PSI - TOLKSL
               MCPDEL  = FXVMC(PPDEL,SNODI(I),PNODI(I),I)
               MCMDEL  = FXVMC(PMDEL,SNODI(I),PNODI(I),I)
               ETAI(I) = (MCPDEL - MCMDEL)/TOL2
               DETAI(I)= (MCPDEL - 2.0D0*MC + MCMDEL)/(TOLKSL*TOLKSL)
               DCKRW(I)= (FXVKR(PPDEL,I) - FXVKR(PMDEL,I))/TOL2
            ELSE
               MCOLD   = FXVMC(POLD,SNODI(I),PNODI(I),I)
               DMCNEW  = FXVDMC(PSI,SNODI(I),PNODI(I),I)
               DMCOLD  = FXVDMC(POLD,SNODI(I),PNODI(I),I)
               ETAI(I) = (MC - MCOLD)/DP
               DETAI(I)= (DMCNEW-DMCOLD)/DP
               DCKRW(I)= (CKRW(I) - FXVKR(POLD,I))/DP
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
               TOL2    = 2.0D0*TOLKSL
               PPDEL   = PSI + TOLKSL
               PMDEL   = PSI - TOLKSL
               SEPDEL  = FHUSE(PPDEL)
               SEMDEL  = FHUSE(PMDEL)
               DSE     = (SEPDEL - SEMDEL)/TOL2
               DDSE    = (SEPDEL - 2.0D0*SE + SEMDEL)/(TOLKSL*TOLKSL)
               DCKRW(I)= (FHUKR2(PPDEL,SEPDEL) - 
     1                   FHUKR2(PMDEL,SEMDEL))/TOL2
            ELSE
               SEOLD   = FHUSE(POLD)
               DSE     = (SE - SEOLD)/DP
               DDSE    = (FHUDSE(PSI) - FHUDSE(POLD))/DP
               DCKRW(I)= (CKRW(I) - FHUKR2(POLD,SEOLD))/DP
            END IF
            DSWDP      = HUSWR1*DSE
            ETAI(I)    = SW(I)*SNODI(I) + PNODI(I)*DSWDP
            DETAI(I)   = DSWDP*SNODI(I) + PNODI(I)*HUSWR1*DDSE
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
               TOL2    = 2.0D0*TOLKSL
               PPDEL   = PSI + TOLKSL
               PMDEL   = PSI - TOLKSL
               SEPDEL  = FHUSE(PPDEL)
               SEMDEL  = FHUSE(PMDEL)
               DSE     = (SEPDEL - SEMDEL)/TOL2
               DDSE    = (SEPDEL - 2.0D0*SE + SEMDEL)/(TOLKSL*TOLKSL)
               DCKRW(I)= (FHUKR3(PPDEL,SEPDEL) - 
     1                   FHUKR3(PMDEL,SEMDEL))/TOL2
            ELSE
               SEOLD   = FHUSE(POLD)
               DSE     = (SE - SEOLD)/DP
               DDSE    = (FHUDSE(PSI) - FHUDSE(POLD))/DP
               DCKRW(I)= (CKRW(I) - FHUKR3(POLD,SEOLD))/DP
            END IF
            DSWDP      = HUSWR1*DSE
            ETAI(I)    = SW(I)*SNODI(I) + PNODI(I)*DSWDP
            DETAI(I)   = DSWDP*SNODI(I) + PNODI(I)*HUSWR1*DDSE
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
               TOL2    = 2.0D0*TOLKSL
               PPDEL   = PSI + TOLKSL
               PMDEL   = PSI - TOLKSL
               SEPDEL  = FBCSE(PPDEL)
               SEMDEL  = FBCSE(PMDEL)
               DSE     = (SEPDEL - SEMDEL)/TOL2
               DDSE    = (SEPDEL - 2.0D0*SE + SEMDEL)/(TOLKSL*TOLKSL)
               DCKRW(I)= (FBCKR(PPDEL) - FBCKR(PMDEL))/TOL2
            ELSE
               SEOLD   = FBCSE(POLD)
               DSE     = (SE - SEOLD)/DP
               DDSE    = (FBCDSE(PSI) - FBCDSE(POLD))/DP
               DCKRW(I)= (CKRW(I) - FBCKR(POLD))/DP
            END IF
            DSWDP      = BCPORM(I)*DSE
            ETAI(I)    = SW(I)*SNODI(I) + PNODI(I)*DSWDP
            DETAI(I)   = DSWDP*SNODI(I) + PNODI(I)*BCPORM(I)*DDSE
         END DO 
      END IF
C
      RETURN
      END
