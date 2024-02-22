C
C**************************  CHPIC4 ************************************
C
C  calculate soil moisture characteristics needed for Picard scheme
C  KSLOPE=4 : localized tangent slope differentiation of moisture 
C  curves.
C  Note that DSETAN contains tangent slope values at each node only for 
C  the case IVGHU=1; for the other IVGHU cases the tangent slope is 
C  constant for all nodes and we use DSETAN(1).
C
C***********************************************************************
C
      SUBROUTINE CHPIC4(N,PTNEW,PNEW,PTIMEP,SWNEW,SWTIMEP,
     1                 SNODI,PNODI,SW,CKRW,ETAI,ET1,ET2,IVGHU)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I
      INTEGER  N,IVGHU
      REAL*8   FVGKR,FVGSE,FVGDSE,FXVKR,FXVMC,FXVDMC
      REAL*8   FHUKR2,FHUKR3,FHUSE,FHUDSE,FBCKR,FBCSE,FBCDSE
      REAL*8   PSI,SE
      REAL*8   PTNEW(*),SNODI(*),PNODI(*),SW(*),CKRW(*),ETAI(*)
      REAL*8   SWNEW(*),SWTIMEP(*),PNEW(*),PTIMEP(*),ET1(*),ET2(*)
      INCLUDE 'SOILCHAR.H'
C
      IF (IVGHU .EQ. 0) THEN
         DO I=1,N
            SWNEW(I)= VGPNOT(I)*FVGSE(PNEW(I),I)+VGRMC(I)/PNODI(I)
            SWTIMEP(I)=VGPNOT(I)*FVGSE(PTIMEP(I),I)+VGRMC(I)/PNODI(I)
         END DO
         DO I=1,N
            PSI       = PTNEW(I)
            SE        = FVGSE(PSI,I)
            SW(I)     = VGPNOT(I)*SE + VGRMC(I)/PNODI(I)
            CKRW(I)   = FVGKR(PSI,SE,I)
            IF (PSI .GE. PSEL  .AND.  PSI .LE. PSER) THEN
              ETAI(I)= SW(I)*SNODI(I) + PNODI(I)*VGPNOT(I)*DSETAN(1)
              ET1(I)=  SW(I)*SNODI(I)
              ET2(I)=  VGPNOT(I)*DSETAN(1)
            ELSE
              ETAI(I)= SW(I)*SNODI(I) + PNODI(I)*VGPNOT(I)*FVGDSE(PSI,I)
              ET1(I)= SW(I)*SNODI(I)
              ET2(I)= VGPNOT(I)*FVGDSE(PSI,I)
            END IF
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
            ET1(I)    = SW(I)*SNODI(I)
            IF (PSI .GE. PSEL  .AND.  PSI .LE. PSER) THEN
               ETAI(I)= DSETAN(I)
               ET2(I) = (ETAI(I)-SW(I)*SNODI(I))/PNODI(I)
            ELSE
               ETAI(I)= FXVDMC(PSI,SNODI(I),PNODI(I),I)
               ET2(I) = (ETAI(I)-SW(I)*SNODI(I))/PNODI(I)
            END IF
         END DO
      ELSE IF (IVGHU .EQ. 2) THEN
       DO I=1,N
          SWNEW(I)=HUSWR1*FHUSE(PNEW(I))+HUSWR
          SWTIMEP(I)=HUSWR1*FHUSE(PTIMEP(I))+HUSWR
       END DO
         DO I=1,N
            PSI       = PTNEW(I)
            SE        = FHUSE(PSI)
            SW(I)     = HUSWR1*SE + HUSWR
            CKRW(I)   = FHUKR2(PSI,SE)
            ET1(I)    = SW(I)*SNODI(I)
            IF (PSI .GE. PSEL  .AND.  PSI .LE. PSER) THEN
               ETAI(I)= SW(I)*SNODI(I) + PNODI(I)*HUSWR1*DSETAN(1)
               ET2(I) = (ETAI(I)-SW(I)*SNODI(I))/PNODI(I)
            ELSE
               ETAI(I)= SW(I)*SNODI(I) + PNODI(I)*HUSWR1*FHUDSE(PSI)
               ET2(I) = (ETAI(I)-SW(I)*SNODI(I))/PNODI(I)
            END IF
         END DO
      ELSE IF (IVGHU .EQ. 3) THEN
       DO I=1,N
          SWNEW(I)=HUSWR1*FHUSE(PNEW(I))+HUSWR
          SWTIMEP(I)=HUSWR1*FHUSE(PTIMEP(I))+HUSWR
       END DO
         DO I=1,N
            PSI       = PTNEW(I)
            SE        = FHUSE(PSI)
            SW(I)     = HUSWR1*SE + HUSWR
            CKRW(I)   = FHUKR3(PSI,SE)
            ET1(I)    = SW(I)*SNODI(I)
            IF (PSI .GE. PSEL  .AND.  PSI .LE. PSER) THEN
               ETAI(I)= SW(I)*SNODI(I) + PNODI(I)*HUSWR1*DSETAN(1)
               ET2(I) = (ETAI(I)-SW(I)*SNODI(I))/PNODI(I)
            ELSE
               ETAI(I)= SW(I)*SNODI(I) + PNODI(I)*HUSWR1*FHUDSE(PSI)
               ET2(I) = (ETAI(I)-SW(I)*SNODI(I))/PNODI(I)
            END IF
         END DO
      ELSE
       DO I=1,N
          SWNEW(I)=BCPORM(I)*FBCSE(PNEW(I))+BCRMC/PNODI(I)
          SWTIMEP(I)=BCPORM(I)*FBCSE(PTIMEP(I))+BCRMC/PNODI(I)
       END DO
         DO I=1,N
            PSI       = PTNEW(I)
            SE        = FBCSE(PSI)
            SW(I)     = BCPORM(I)*SE + BCRMC/PNODI(I)
            CKRW(I)   = FBCKR(PSI)
            ET1(I)    = SW(I)*SNODI(I)
            IF (PSI .GE. PSEL  .AND.  PSI .LE. PSER) THEN
               ETAI(I)= SW(I)*SNODI(I) + PNODI(I)*BCPORM(I)*DSETAN(1)
               ET2(I) =(ETAI(I)-SW(I)*SNODI(I))/PNODI(I)
            ELSE
               ETAI(I)= SW(I)*SNODI(I) + PNODI(I)*BCPORM(I)*FBCDSE(PSI)
               ET2(I) =(ETAI(I)-SW(I)*SNODI(I))/PNODI(I)
            END IF
         END DO
      END IF
C
      RETURN
      END
