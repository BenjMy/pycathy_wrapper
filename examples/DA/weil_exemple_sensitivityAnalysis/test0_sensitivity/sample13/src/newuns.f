C
C**************************  NEWUNS ************************************
C
C  calculate characteristic curves for unsaturated zone,
C  Newton case
C
C***********************************************************************
C
      SUBROUTINE NEWUNS(NLKP,N,NT,KSLOPE,IVGHU,NTRI,TP,TETRA,
     1                  PTNEW,PTOLD,SNODI,PNODI,
     2                  SW,CKRW,ETAI,DCKRW,DETAI,
     3                  SWE,CKRWE,ETAE,DCKRWE,DETAIE)
C
      IMPLICIT NONE
      INTEGER  NLKP,N,NT,KSLOPE,IVGHU,NTRI
      INTEGER  TP(*),TETRA(5,*)
      REAL*8   PTNEW(*),PTOLD(*),SNODI(*),PNODI(*)
      REAL*8   SW(*),CKRW(*),ETAI(*),DCKRW(*),DETAI(*)
      REAL*8   SWE(*),CKRWE(*),ETAE(*),DCKRWE(*),DETAIE(*)
C
      IF (IVGHU .EQ. -1) THEN
         CALL MOISTABNEW(NLKP,NT,NTRI,TETRA,PTNEW,PNODI,SNODI,
     1                   SWE,CKRWE,ETAE,DCKRWE,DETAIE)
         CALL ELTNOD(N,NT,TP,TETRA,DETAIE,DETAI)
         CALL ELTNOD(N,NT,TP,TETRA,DCKRWE,DCKRW)
         CALL ELTNOD(N,NT,TP,TETRA,ETAE,ETAI)
         CALL ELTNOD(N,NT,TP,TETRA,CKRWE,CKRW)
         CALL ELTNOD(N,NT,TP,TETRA,SWE,SW)
      ELSE
         IF (KSLOPE .EQ. 0) THEN
            CALL CHNEW0(N,PTNEW,SNODI,PNODI,SW,CKRW,ETAI,
     1                  DCKRW,DETAI,IVGHU)
         ELSE IF (KSLOPE .EQ. 1) THEN
            CALL CHNEW1(N,PTNEW,PTOLD,SNODI,PNODI,SW,CKRW,ETAI,
     1                  DCKRW,DETAI,IVGHU)
         ELSE IF (KSLOPE .EQ. 2) THEN
            CALL CHNEW2(N,PTNEW,PTOLD,SNODI,PNODI,SW,CKRW,ETAI,
     1                  DCKRW,DETAI,IVGHU)
         ELSE IF (KSLOPE .EQ. 3) THEN
            CALL CHNEW3(N,PTNEW,PTOLD,SNODI,PNODI,SW,CKRW,ETAI,
     1                  DCKRW,DETAI,IVGHU)
         ELSE
            CALL CHNEW4(N,PTNEW,SNODI,PNODI,SW,CKRW,ETAI,DCKRW,DETAI,
     1                  IVGHU)
         END IF
         CALL NODELT(NT,TETRA,CKRW,CKRWE)
         CALL NODELT(NT,TETRA,ETAI,ETAE)
      END IF
C
      RETURN
      END
