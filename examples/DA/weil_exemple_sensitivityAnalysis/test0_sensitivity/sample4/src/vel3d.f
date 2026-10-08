C
C**************************  VEL3D  ************************************
C
C  subroutine per il calcolo delle velocita' di Darcy lungo l'asse
C  X - Y - Z nel baricentro di ogni tetraedro
C
C  UU = velocita' di Darcy uniforme nel tetraedro lungo X
C  VV = velocita' di Darcy uniforme nel tetraedro lungo Y
C  WW = velocita' di Darcy uniforme nel tetraedro lungo Z
C
C***********************************************************************
C
      SUBROUTINE VEL3D(NT,TETRA,NTRI,PNEW,PERMX,PERMY,PERMZ,
     1                 UU,VV,WW,BI,CI,DI,CKRWE,VOLUR,IVOL)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I,J,KK,IR,ISTR,MTYPE
      INTEGER  NT,NTRI
      INTEGER  TETRA(5,*),IVOL(*)
      REAL*8   BB,CC,DD,XYZ
      REAL*8   PERMX(MAXSTR,*),PERMY(MAXSTR,*),PERMZ(MAXSTR,*)
      REAL*8   PNEW(*),CKRWE(*)
      REAL*8   UU(*),VV(*),WW(*),BI(4,*),CI(4,*),DI(4,*),VOLUR(*)
C
C  inizia il ciclo sui tetraedri
C
      DO KK=1,NT
C
C  ISTR = numero dello strato a cui appartiene il tetraedro KK
C
         ISTR=1+KK/(NTRI*3)
         IR=MOD(KK,NTRI*3)
         IF(IR.EQ.0) ISTR=ISTR-1
         MTYPE=TETRA(5,KK)
C
C  calcola la velocita' nel tetraedro
C
         BB=0.0D0
         CC=0.0D0
         DD=0.0D0
         XYZ=-CKRWE(KK)*VOLUR(KK)*IVOL(KK)
         DO J=1,4
            I=TETRA(J,KK)
            BB=BB+PNEW(I)*BI(J,KK)
            CC=CC+PNEW(I)*CI(J,KK)
            DD=DD+PNEW(I)*DI(J,KK)
         END DO
         UU(KK)=BB*XYZ*PERMX(ISTR,MTYPE)
         VV(KK)=CC*XYZ*PERMY(ISTR,MTYPE)
         WW(KK)=(DD*XYZ - CKRWE(KK))*PERMZ(ISTR,MTYPE)
      END DO
C
      RETURN
      END
