C
C**************************  CHPIC2 ************************************
C
C  calculate soil moisture characteristics needed for Picard scheme
C  KSLOPE=2 : chord slope and centered difference formulas for 
C             differentiation of moisture curves
C
C***********************************************************************
C
      SUBROUTINE CHPIC2(N,PTNEW,PTOLD,PNEW,PTIMEP,SWNEW,SWTIMEP,
     1                  SNODI,PNODI,SW,CKRW,ETAI,ET1,ET2,IVGHU)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I
      INTEGER  N,IVGHU
      REAL*8   FVGKR,FVGSE,FXVKR,FXVMC
      REAL*8   FHUKR2,FHUKR3,FHUSE,FBCKR,FBCSE
      REAL*8   PSI,POLD,DP,PPDEL,PMDEL,SE,DSE,MC,MCOLD,MCPDEL,MCMDEL
      REAL*8   PTNEW(*),PTOLD(*),SNODI(*),PNODI(*),SW(*),CKRW(*),ETAI(*)
      REAL*8   ET1(*),ET2(*),SWNEW(*),SWTIMEP(*),PNEW(*),PTIMEP(*)
      INCLUDE 'SOILCHAR.H'
C
      IF (IVGHU .EQ. 0) THEN
         DO I=1,N
            SWNEW(I)=VGPNOT(I)*FVGSE(PNEW(I),I)+VGRMC(I)/PNODI(I)
            SWTIMEP(I)=VGPNOT(I)*FVGSE(PTIMEP(I),I)+VGRMC(I)/PNODI(I)
         END DO
         DO I=1,N
            PSI      = PTNEW(I)
            SE       = FVGSE(PSI,I)
            SW(I)    = VGPNOT(I)*SE + VGRMC(I)/PNODI(I)
            CKRW(I)  = FVGKR(PSI,SE,I)
            POLD     = PTOLD(I)
            DP       = PSI - POLD
            IF (DABS(DP) .LT. TOLKSL) THEN
               PPDEL = PSI + TOLKSL
               PMDEL = PSI - TOLKSL
               DSE   = (FVGSE(PPDEL,I) - FVGSE(PMDEL,I))/(2.0D0*TOLKSL)
            ELSE
               DSE   = (SE - FVGSE(POLD,I))/DP
            END IF
            ETAI(I)  = SW(I)*SNODI(I) + PNODI(I)*VGPNOT(I)*DSE
            ET1(I)   = SW(I)*SNODI(I)
            ET2(I)   = VGPNOT(I)*DSE
         END DO
      ELSE IF (IVGHU .EQ. 1) THEN
       DO I=1,N
          SWNEW(I)=FXVMC(PNEW(I),SNODI(I),PNODI(I),I)/PNODI(I)
          SWTIMEP(I)=FXVMC(PTIMEP(I),SNODI(I),PNODI(I),I)/PNODI(I)
       END DO 
         DO I=1,N
            PSI       = PTNEW(I)
            SW(I)     = FXVMC(PSI,SNODI(I),PNODI(I),I)/PNODI(I)
            CKRW(I)   = FXVKR(PSI,I)
            POLD      = PTOLD(I)
            DP        = PSI - POLD
            IF (DABS(DP) .LT. TOLKSL) THEN
               PPDEL  = PSI + TOLKSL
               PMDEL  = PSI - TOLKSL
               MCPDEL = FXVMC(PPDEL,SNODI(I),PNODI(I),I)
               MCMDEL = FXVMC(PMDEL,SNODI(I),PNODI(I),I)
               ETAI(I)= (MCPDEL - MCMDEL)/(2.0D0*TOLKSL)
               ET1(I) = SW(I)*SNODI(I) 
               ET2(I) = (ETAI(I)-SNODI(I)*SW(I))/PNODI(I)
            ELSE
               MC     = FXVMC(PSI,SNODI(I),PNODI(I),I)
               MCOLD  = FXVMC(POLD,SNODI(I),PNODI(I),I)
               ETAI(I)= (MC - MCOLD)/DP
               ET1(I) = SW(I)*SNODI(I)
               ET2(I) = (ETAI(I)-SW(I)*SNODI(I))/PNODI(I)
            END IF
         END DO
      ELSE IF (IVGHU .EQ. 2) THEN
       DO I=1,N
          SWNEW(I)=HUSWR1*FHUSE(PNEW(I))+HUSWR
          SWTIMEP(I)=HUSWR1*FHUSE(PTIMEP(I))+HUSWR
       END DO
         DO I=1,N
            PSI     = PTNEW(I)
            SE      = FHUSE(PSI)
            SW(I)   = HUSWR1*SE + HUSWR
            CKRW(I) = FHUKR2(PSI,SE)
            POLD    = PTOLD(I)
            DP      = PSI - POLD
            IF (DABS(DP) .LT. TOLKSL) THEN
               PPDEL= PSI + TOLKSL
               PMDEL= PSI - TOLKSL
               DSE  = (FHUSE(PPDEL) - FHUSE(PMDEL))/(2.0D0*TOLKSL)
            ELSE
               DSE  = (SE - FHUSE(POLD))/DP
            END IF
            ETAI(I) = SW(I)*SNODI(I) + PNODI(I)*HUSWR1*DSE
            ET1(I) = SW(I)*SNODI(I)
            ET2(I) = (ETAI(I)-SW(I)*SNODI(I))/PNODI(I) 
         END DO
      ELSE IF (IVGHU .EQ. 3) THEN
       DO I=1,N
          SWNEW(I)=HUSWR1*FHUSE(PNEW(I))+HUSWR
          SWTIMEP(I)=HUSWR1*FHUSE(PTIMEP(I))+HUSWR
       END DO
         DO I=1,N
            PSI     = PTNEW(I)
            SE      = FHUSE(PSI)
            SW(I)   = HUSWR1*SE + HUSWR
            CKRW(I) = FHUKR3(PSI,SE)
            POLD    = PTOLD(I)
            DP      = PSI - POLD
            IF (DABS(DP) .LT. TOLKSL) THEN
               PPDEL= PSI + TOLKSL
               PMDEL= PSI - TOLKSL
               DSE  = (FHUSE(PPDEL) - FHUSE(PMDEL))/(2.0D0*TOLKSL)
            ELSE
               DSE  = (SE - FHUSE(POLD))/DP
            END IF
            ETAI(I) = SW(I)*SNODI(I) + PNODI(I)*HUSWR1*DSE
            ET1(I)  = SW(I)*SNODI(I)
            ET2(I)  = (ETAI(I)-SW(I)*SNODI(I))/PNODI(I)
         END DO
      ELSE
       DO I=1,N
          SWNEW(I)=BCPORM(I)*FBCSE(PNEW(I))+BCRMC/PNODI(I)
          SWTIMEP(I)=BCPORM(I)*FBCSE(PTIMEP(I))+BCRMC/PNODI(I)
       END DO
         DO I=1,N
            PSI      = PTNEW(I)
            SE       = FBCSE(PSI)
            SW(I)    = BCPORM(I)*SE + BCRMC/PNODI(I)
            CKRW(I)  = FBCKR(PSI)
            POLD     = PTOLD(I)
            DP       = PSI - POLD
            IF (DABS(DP) .LT. TOLKSL) THEN
               PPDEL = PSI + TOLKSL
               PMDEL = PSI - TOLKSL
               DSE   = (FBCSE(PPDEL) - FBCSE(PMDEL))/(2.0D0*TOLKSL)
            ELSE
               DSE   = (SE - FBCSE(POLD))/DP
            END IF
            ETAI(I)  = SW(I)*SNODI(I) + PNODI(I)*BCPORM(I)*DSE
            ET1(I)   = SW(I)*SNODI(I)
            ET2(I)   = (ETAI(I)-SW(I)*SNODI(I))/PNODI(I)
         END DO
      END IF
C
      RETURN
      END
