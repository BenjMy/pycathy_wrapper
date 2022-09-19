C
C**************************  CHPIC0 ************************************
C
C  calculate soil moisture characteristics needed for Picard scheme
C  KSLOPE=0 : analytical differentiation of moisture curves
C
C***********************************************************************
C
      SUBROUTINE CHPIC0(N,PTNEW,PNEW,PTIMEP,SWNEW,SWTIMEP,
     1                  SNODI,PNODI,SW,CKRW,ETAI,ET1,ET2,IVGHU)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I
      INTEGER  N,IVGHU
      REAL*8   FVGKR,FVGSE,FVGDSE,FXVKR,FXVMC,FXVDMC
      REAL*8   FHUKR2,FHUKR3,FHUSE,FHUDSE,FBCKR,FBCSE,FBCDSE
      REAL*8   PSI,SE
      REAL*8   PTNEW(*),SNODI(*),PNODI(*),SW(*),CKRW(*),ETAI(*)
      REAL*8   ET1(*),ET2(*),PNEW(*),PTIMEP(*),SWNEW(*),SWTIMEP(*)
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
            ETAI(I)  = SW(I)*SNODI(I) + PNODI(I)*VGPNOT(I)*FVGDSE(PSI,I)
            ET1(I)   = SW(I)*SNODI(I)
            ET2(I)   = VGPNOT(I)*FVGDSE(PSI,I)
            CKRW(I)  = FVGKR(PSI,SE,I)
         END DO
      ELSE IF (IVGHU .EQ. 1) THEN
       DO I=1,N
          SWNEW(I)=FXVMC(PNEW(I),SNODI(I),PNODI(I),I)/PNODI(I)
          SWTIMEP(I)=FXVMC(PTIMEP(I),SNODI(I),PNODI(I),I)/PNODI(I)
       END DO
         DO I=1,N
            PSI      = PTNEW(I)
            SW(I)    = FXVMC(PSI,SNODI(I),PNODI(I),I)/PNODI(I)
            ETAI(I)  = FXVDMC(PSI,SNODI(I),PNODI(I),I)
            ET1(I)   = SW(I)*SNODI(I)
            ET2(I)   = (ETAI(I)-SW(I)*SNODI(I))/PNODI(I)
            CKRW(I)  = FXVKR(PSI,I)
         END DO
      ELSE IF (IVGHU .EQ. 2) THEN
       DO I=1,N
          SWNEW(I)=HUSWR1*FHUSE(PNEW(I))+HUSWR
          SWTIMEP(I)=HUSWR1*FHUSE(PTIMEP(I))+HUSWR
       END DO
         DO I=1,N
            PSI      = PTNEW(I)
            SE       = FHUSE(PSI)
            SW(I)    = HUSWR1*SE + HUSWR
            ETAI(I)  = SW(I)*SNODI(I) + PNODI(I)*HUSWR1*FHUDSE(PSI)
            ET1(I)  =  SW(I)*SNODI(I)
            ET2(I)   = (ETAI(I)-SW(I)*SNODI(I))/PNODI(I) 
            CKRW(I)  = FHUKR2(PSI,SE)
         END DO
      ELSE IF (IVGHU .EQ. 3) THEN
       DO I=1,N
          SWNEW(I)=HUSWR1*FHUSE(PNEW(I))+HUSWR
          SWTIMEP(I)=HUSWR1*FHUSE(PTIMEP(I))+HUSWR
       END DO
         DO I=1,N
            PSI      = PTNEW(I)
            SE       = FHUSE(PSI)
            SW(I)    = HUSWR1*SE + HUSWR
            ETAI(I)  = SW(I)*SNODI(I) + PNODI(I)*HUSWR1*FHUDSE(PSI)
            ET1(I)  = SW(I)*SNODI(I)
            ET2(I)  = (ETAI(I)-SW(I)*SNODI(I))/PNODI(I)
            CKRW(I)  = FHUKR3(PSI,SE)
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
            ETAI(I)  = SW(I)*SNODI(I) + PNODI(I)*BCPORM(I)*FBCDSE(PSI)
            ET1(I)   = SW(I)*SNODI(I)
            ET2(I)   = (ETAI(I)-SW(I)*SNODI(I))/PNODI(I)
            CKRW(I)  = FBCKR(PSI)
         END DO
      END IF
C
      RETURN
      END
