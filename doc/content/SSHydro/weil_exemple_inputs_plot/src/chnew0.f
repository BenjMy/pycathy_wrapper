C
C**************************  CHNEW0 ************************************
C
C  calculate soil moisture characteristics needed for Newton scheme
C  KSLOPE=0 : analytical differentiation of moisture curves
C
C***********************************************************************
C
      SUBROUTINE CHNEW0(N,PTNEW,SNODI,PNODI,SW,CKRW,ETAI,
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
      REAL*8   PSI,SE,DSEDP,DSWDP
      REAL*8   PTNEW(*),SNODI(*),PNODI(*)
      REAL*8   SW(*),CKRW(*),ETAI(*),DCKRW(*),DETAI(*)
      INCLUDE 'SOILCHAR.H'
C
      IF (IVGHU .EQ. 0) THEN
         DO I=1,N
            PSI      = PTNEW(I)
            SE       = FVGSE(PSI,I)
            DSEDP    = FVGDSE(PSI,I)
            DSWDP    = VGPNOT(I)*DSEDP
            SW(I)    = VGPNOT(I)*SE + VGRMC(I)/PNODI(I)
            ETAI(I)  = SW(I)*SNODI(I) + PNODI(I)*DSWDP
            DETAI(I) = DSWDP*SNODI(I) + PNODI(I)*VGPNOT(I)*FVGDDS(PSI,I)
            CKRW(I)  = FVGKR(PSI,SE,I)
            DCKRW(I) = FVGDKR(PSI,I)
         END DO 
      ELSE IF (IVGHU .EQ. 1) THEN
         DO I=1,N
            PSI      = PTNEW(I)
            SW(I)    = FXVMC(PSI,SNODI(I),PNODI(I),I)/PNODI(I)
            ETAI(I)  = FXVDMC(PSI,SNODI(I),PNODI(I),I)
            DETAI(I) = FXVDDM(PSI,PNODI(I),I)
            CKRW(I)  = FXVKR(PSI,I)
            DCKRW(I) = FXVDKR(PSI,I)
         END DO
      ELSE IF (IVGHU .EQ. 2) THEN
         DO I=1,N
            PSI      = PTNEW(I)
            SE       = FHUSE(PSI)
            DSEDP    = FHUDSE(PSI)
            DSWDP    = HUSWR1*DSEDP
            SW(I)    = HUSWR1*SE + HUSWR
            ETAI(I)  = SW(I)*SNODI(I) + PNODI(I)*DSWDP
            DETAI(I) = DSWDP*SNODI(I) + PNODI(I)*HUSWR1*FHUDDS(PSI)
            CKRW(I)  = FHUKR2(PSI,SE)
            DCKRW(I) = FHUDK2(PSI,SE,DSEDP)
         END DO
      ELSE IF (IVGHU .EQ. 3) THEN
         DO I=1,N
            PSI      = PTNEW(I)
            SE       = FHUSE(PSI)
            DSEDP    = FHUDSE(PSI)
            DSWDP    = HUSWR1*DSEDP
            SW(I)    = HUSWR1*SE + HUSWR
            ETAI(I)  = SW(I)*SNODI(I) + PNODI(I)*DSWDP
            DETAI(I) = DSWDP*SNODI(I) + PNODI(I)*HUSWR1*FHUDDS(PSI)
            CKRW(I)  = FHUKR3(PSI,SE)
            DCKRW(I) = FHUDK3(PSI,SE,DSEDP,CKRW(I))
         END DO
      ELSE
         DO I=1,N
            PSI      = PTNEW(I)
            SE       = FBCSE(PSI)
            DSEDP    = FBCDSE(PSI)
            DSWDP    = BCPORM(I)*DSEDP
            SW(I)    = BCPORM(I)*SE + BCRMC/PNODI(I)
            ETAI(I)  = SW(I)*SNODI(I) + PNODI(I)*DSWDP
            DETAI(I) = DSWDP*SNODI(I) + PNODI(I)*BCPORM(I)*FBCDDS(PSI)
            CKRW(I)  = FBCKR(PSI)
            DCKRW(I) = FBCDKR(PSI)
         END DO 
      END IF
C
      RETURN
      END
