C
C**************************  CHPARM ************************************
C
C  calculate miscellaneous van Genuchten, Huyakorn, and
C  Brooks-Corey moisture curve parameters, including VGPNOT for 
C  van Genuchten and extended van Genuchten cases and BCPORM
C  for Brooks-Corey case
C
C***********************************************************************
C
      SUBROUTINE CHPARM(N,SNODI,PNODI,IVGHU)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I,K,IK
      INTEGER  N,IVGHU
      REAL*8   SS,TSR,VGM1,PNOT,V1,V2,V3,V4,V5,V6,V7,V8
      REAL*8   SNODI(*),PNODI(*)
      INCLUDE 'SOILCHAR.H'
      INCLUDE 'IOUNITS.H'
C
      IF (IVGHU .EQ. 0 .OR. IVGHU .EQ. 1) THEN
         DO I=1,N
            VGM(I)=(VGN(I)-1.0D0)/VGN(I)
            VGN1(I)=VGN(I)-1.0D0
            VGNR(I)=1.0D0/VGN(I)
            VGMM1(I)=VGM(I)-1.0D0
            VGPSN(I)=DABS(VGPSAT(I))**VGN(I)
            VGMR(I)=1.0D0/VGM(I)
            VGM52(I)=2.5D0*VGM(I)
         END DO
         IF (IVGHU .EQ. 0) THEN
            DO I=1,N
               VGPNOT(I)=(PNODI(I)-VGRMC(I))/PNODI(I)
            END DO
         ELSE
C
C  bisection algorithm to compute VGPNOT for extended van Genuchten case
C
            IK=500
            V1=1.0D-14
            DO I=1,N
               VGM1=VGM(I)+1.0D0 
               SS=SNODI(I)
               TSR=PNODI(I)-VGRMC(I)
               V2=VGPSAT(I)*(VGM(I)**(1.0D0/VGN(I)))
               V3=0.0D0
               V4=-VGM(I)*VGN(I)*TSR*(VGM(I)**VGM(I))/
     1             (VGPSAT(I)*(VGM1**VGM1))
               IF (SS .GE. V4) THEN
                  WRITE(IOUT2,900) I,SS,V4
                  CALL CLOSIO
                  STOP
               END IF
               K=0
   30          K=K+1
                  PNOT=V2 + (V3-V2)/2.0D0
                  V5=(PNOT/VGPSAT(I))**VGN(I)
                  V6=(1.0D0+V5)**VGM1
                  V7=((DABS(PNOT)**VGN1(I))/V6) - (SS*VGPSN(I)/
     1               (VGN1(I)*TSR))
                  IF (V7 .EQ. 0.0D0 .OR. (V3-V2)/2.0D0 .LT. V1) GO TO 50
                  V5=(V2/VGPSAT(I))**VGN(I)
                  V6=(1.0D0+V5)**VGM1
                  V8=((DABS(V2)**VGN1(I))/V6) - (SS*VGPSN(I)/
     1               (VGN1(I)*TSR))
                  IF (V8*V7 .GT. 0.0D0) THEN
                     V2=PNOT
                  ELSE
                     V3=PNOT
                  END IF
                  IF (K .LT. IK) GO TO 30
                  WRITE(IOUT2,910) IK
                  CALL CLOSIO
                  STOP
   50          VGPNOT(I)=PNOT
            END DO
         END IF
      ELSE IF (IVGHU .EQ. 2 .OR. IVGHU .EQ. 3) THEN
         HUSWR1=1.0D0-HUSWR
         HUALB=HUALFA**int(HUBETA)
         HUGAM1=HUGAMA+1.0D0
         HUGB=HUGAMA*HUBETA
         HUGB1=1.0D0+HUGB
         HU1BET=1.0D0-HUBETA
         HUGAM2=HUGAMA+2.0D0
         IF (IVGHU .EQ. 2) THEN
            HUN1=HUN-1.0D0
         ELSE
            HU2A=2.0D0*HUA
            HUB2A=HUB-HU2A
            HUAB=HUA-HUB
            HULN10=DLOG(10.0D0)
         END IF
      ELSE IF (IVGHU .EQ. 4) THEN
         DO I=1,N
            BCPORM(I)=(PNODI(I)-BCRMC)/PNODI(I)
         END DO
         BCB1=BCBETA+1.0D0
         BCB2=BCBETA+2.0D0
         BCBB1P=BCBETA*BCB1/(BCPSAT*BCPSAT)
         BCBPS=BCBETA/DABS(BCPSAT)
         BC23B=2.0D0 + (3.0D0*BCBETA)
         BC23BP=BC23B/DABS(BCPSAT)
         BC33B=3.0D0 + (3.0D0*BCBETA)
      END IF
C
      RETURN
  900 FORMAT(/,'  ****** ERROR IN SNODI INPUT: NODE ',I5,
     1         ';   SNODI = ',1PE13.5,
     2       /,'  ****** SS MUST BE SMALLER THAN DMCMAX = ',1PE13.5,
     3       /,'  ****** PROGRAM ABORTING')
  910 FORMAT(/,'  ****** MAX ITERATIONS = ',I5,' EXCEEDED IN BISECTION',
     1         ' ALGORITHM',
     2       /,'  ****** PROGRAM ABORTING')
      END
