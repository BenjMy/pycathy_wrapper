C
C**************************  PICUNS ************************************
C
C  calculate characteristic curves for unsaturated zone,
C  Picard case
C
C***********************************************************************
C
      SUBROUTINE PICUNS(NLKP,N,NT,KSLOPE,IVGHU,NTRI,TP,TETRA,
     1                  PTNEW,PTOLD,PNEW,PTIMEP,SWNEW,SWTIMEP,
     2                  SNODI,PNODI,SW,CKRW,ETAI,ET1,ET2,ET1E,  
     2                  SWE,CKRWE,ETAE,PEL)
C
      IMPLICIT NONE
      INTEGER  NLKP,N,NT,KSLOPE,IVGHU,NTRI
      INTEGER  TP(*),TETRA(5,*)
      REAL*8   PTNEW(*),PTOLD(*),SNODI(*),PNODI(*)
      REAL*8   PNEW(*),PTIMEP(*),SWNEW(*),SWTIMEP(*)
      REAL*8   SW(*),CKRW(*),ETAI(*),SWE(*),CKRWE(*),ETAE(*)
      REAL*8   ET1(*),ET2(*),ET1E(*),PEL(*)
C
      IF (IVGHU .EQ. -1) THEN
         CALL MOISTABPIC(NLKP,NT,NTRI,TETRA,PTNEW,PNODI,SNODI,
     1                   SWE,CKRWE,ETAE)
         CALL ELTNOD(N,NT,TP,TETRA,ETAE,ETAI)
         CALL ELTNOD(N,NT,TP,TETRA,CKRWE,CKRW)
         CALL ELTNOD(N,NT,TP,TETRA,SWE,SW)
      ELSE
         IF (KSLOPE .EQ. 0) THEN
            CALL CHPIC0(N,PTNEW,PNEW,PTIMEP,SWNEW,SWTIMEP,
     1                  SNODI,PNODI,SW,CKRW,ETAI,ET1,ET2,IVGHU)
         ELSE IF (KSLOPE .EQ. 1) THEN
            CALL CHPIC1(N,PTNEW,PTOLD,PNEW,PTIMEP,SWNEW,SWTIMEP,
     1                 SNODI,PNODI,SW,CKRW,ETAI,ET1,ET2,IVGHU)
         ELSE IF (KSLOPE .EQ. 2) THEN
            CALL CHPIC2(N,PTNEW,PTOLD,PNEW,PTIMEP,SWNEW,SWTIMEP,
     1                 SNODI,PNODI,SW,CKRW,ETAI,ET1,ET2,IVGHU)
         ELSE IF (KSLOPE .EQ. 3) THEN
            CALL CHPIC3(N,PTNEW,PTOLD,PNEW,PTIMEP,SWNEW,SWTIMEP,
     1                  SNODI,PNODI,SW,CKRW,ETAI,ET1,ET2,IVGHU)
         ELSE
            CALL CHPIC4(N,PTNEW,PNEW,PTIMEP,SWNEW,SWTIMEP,
     1                  SNODI,PNODI,SW,CKRW,ETAI,ET1,ET2,IVGHU)
         END IF
         CALL NODELT(NT,TETRA,CKRW,CKRWE)
         CALL NODELT(NT,TETRA,ETAI,ETAE)
         CALL NODELT(NT,TETRA,ET1,ET1E)
         CALL NODELT(NT,TETRA,PNODI,PEL)
      END IF
C
      RETURN
      END
