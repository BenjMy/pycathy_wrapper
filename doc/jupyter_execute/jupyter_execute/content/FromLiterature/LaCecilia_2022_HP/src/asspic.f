C
C**************************  ASSPIC ************************************
C
C  assemble global stiffness and mass matrices from the local
C  contributions of each element: Picard scheme
C
C***********************************************************************
C
      SUBROUTINE ASSPIC(NT,NTRI,TETRA,TETJA,LMASS,COEF1,COEF2,
     1                  COEF4,ET1E,PEL,
     1                  PERMX,PERMY,PERMZ,CKRWE,ETAE,
     2                  VOLU,VOLUR,BI,CI,DI)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  K,L,IEL,ISTR,IR,IND,MTYPE
      INTEGER  NT,NTRI
      INTEGER  TETRA(5,*),TETJA(4,4,*)
      REAL*8   VK,VE,VE1,PVK1,PVK2,PVK3,PVKB,PVKC,PVKD
      REAL*8   CKRWE(*),ETAE(*),VOLU(*),VOLUR(*)
      REAL*8   ET1E(*),PEL(*)
      REAL*8   COEF1(*),COEF2(*),COEF4(*)
      REAL*8   PERMX(MAXSTR,*),PERMY(MAXSTR,*),PERMZ(MAXSTR,*)
      REAL*8   LMASS(4,4),BI(4,*),CI(4,*),DI(4,*)
C
      DO IEL=1,NT
         ISTR=1+IEL/(NTRI*3)
         IR=MOD(IEL,NTRI*3)
         IF (IR .EQ. 0) ISTR=ISTR-1
         MTYPE=TETRA(5,IEL)
         VK=VOLUR(IEL)*CKRWE(IEL)
c CARLOTTA I separate the two coef for the storage
         VE1=VOLU(IEL)*ET1E(IEL)
         VE=VOLU(IEL)*PEL(IEL)
c  CARLOTTA
         PVK1=PERMX(ISTR,MTYPE)*VK
         PVK2=PERMY(ISTR,MTYPE)*VK
         PVK3=PERMZ(ISTR,MTYPE)*VK
         DO K=1,4
            PVKB=PVK1*BI(K,IEL)
            PVKC=PVK2*CI(K,IEL)
            PVKD=PVK3*DI(K,IEL)
            DO L=K,4
               IND=TETJA(K,L,IEL)
               COEF1(IND)=COEF1(IND) + PVKB*BI(L,IEL) + PVKC*CI(L,IEL) +
     1                                 PVKD*DI(L,IEL)
               COEF2(IND)=COEF2(IND) + VE1*LMASS(K,L)
               COEF4(IND)=COEF4(IND) + VE*LMASS(K,L)
            END DO
         END DO
      END DO
C
      RETURN
      END
