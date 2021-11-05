C
C**************************  CHVELO ************************************
C
C  calculate soil moisture characteristics needed for storage and
C  velocity calculations and for output 
C
C***********************************************************************
C
      SUBROUTINE CHVELO(NLKP,N,NT,NTRI,IVGHU,TETRA,TP,
     1                  PNEW,SNODI,PNODI,SW,CKRW,SWE,CKRWE)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I
      INTEGER  NLKP,N,NT,NTRI,IVGHU
      INTEGER  TETRA(5,*),TP(*)
      REAL*8   FVGKR,FVGSE,FXVKR,FXVMC
      REAL*8   FHUKR2,FHUKR3,FHUSE,FBCKR,FBCSE
      REAL*8   PSI,SE
      REAL*8   PNEW(*),SNODI(*),PNODI(*),SW(*),CKRW(*)
      REAL*8   SWE(*),CKRWE(*)
      INCLUDE 'SOILCHAR.H'
C
      IF (IVGHU .EQ. -1) THEN
         CALL MOISTABVEL(NLKP,NT,NTRI,TETRA,PNEW,SWE,CKRWE)
         CALL ELTNOD(N,NT,TP,TETRA,CKRWE,CKRW)
         CALL ELTNOD(N,NT,TP,TETRA,SWE,SW)
      ELSE IF (IVGHU .EQ. 0) THEN
         DO I=1,N
            PSI      = PNEW(I)
            SE       = FVGSE(PSI)
            SW(I)    = VGPNOT(I)*SE + VGRMC/PNODI(I)
            CKRW(I)  = FVGKR(PSI,SE)
         END DO 
      ELSE IF (IVGHU .EQ. 1) THEN
         DO I=1,N
            PSI      = PNEW(I)
            SW(I)    = FXVMC(PSI,SNODI(I),PNODI(I),I)/PNODI(I)
            CKRW(I)  = FXVKR(PSI)
         END DO
      ELSE IF (IVGHU .EQ. 2) THEN
         DO I=1,N
            PSI      = PNEW(I)
            SE       = FHUSE(PSI)
            SW(I)    = HUSWR1*SE + HUSWR
            CKRW(I)  = FHUKR2(PSI,SE)
         END DO
      ELSE IF (IVGHU .EQ. 3) THEN
         DO I=1,N
            PSI      = PNEW(I)
            SE       = FHUSE(PSI)
            SW(I)    = HUSWR1*SE + HUSWR
            CKRW(I)  = FHUKR3(PSI,SE)
         END DO
      ELSE 
         DO I=1,N
            PSI      = PNEW(I)
            SE       = FBCSE(PSI)
            SW(I)    = BCPORM(I)*SE + BCRMC/PNODI(I)
            CKRW(I)  = FBCKR(PSI)
         END DO 
      END IF
C
      RETURN
      END
