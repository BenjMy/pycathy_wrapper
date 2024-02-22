C
C**************************  ASSNEW ************************************
C
C  assemble global stiffness and mass matrices from the local
C  contributions of each element, and also the derivative term
C  components of the Jacobian: Newton scheme
C
C***********************************************************************
C
      SUBROUTINE ASSNEW(NT,NTRI,TETRA,TETJA,LMASS,COEF1,COEF2,
     1                  COEF3,PERMX,PERMY,PERMZ,CKRWE,ETAE,
     2                  PTNEW,PTIMEP,PNEW,DCKRW,DETAI,TETAF,DELTAT,
     3                  VOLU,VOLUR,IVOL,BI,CI,DI)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  K,L,M,IEL,ISTR,IR,KNOD,LNOD,MNOD,IND,INDS,MTYPE
      INTEGER  NT,NTRI
      INTEGER  TETRA(5,*),TETJA(4,4,*),IVOL(*)
      REAL*8   RDT,VE,PVK1,PVK2,PVK3,PVKB,PVKC,PVKD,TSUMTD,SUM1TV
      REAL*8   PV1,PV2,PV3,PVB,PVC,PVD,TPI,TVD,SUM,SUM1
      REAL*8   TETAF,DELTAT
      REAL*8   CKRWE(*),ETAE(*),VOLU(*),VOLUR(*)
      REAL*8   COEF1(*),COEF2(*),COEF3(*)
      REAL*8   PERMX(MAXSTR,*),PERMY(MAXSTR,*),PERMZ(MAXSTR,*)
      REAL*8   PNEW(*),PTNEW(*),PTIMEP(*),DCKRW(*),DETAI(*)
      REAL*8   LMASS(4,4),BI(4,*),CI(4,*),DI(4,*)
C
      RDT=1.0D0/DELTAT
      DO IEL=1,NT
         ISTR=1+IEL/(NTRI*3)
         IR=MOD(IEL,NTRI*3)
         IF (IR .EQ. 0) ISTR=ISTR-1
         MTYPE=TETRA(5,IEL)
         VE=VOLU(IEL)*ETAE(IEL)
         PV1=PERMX(ISTR,MTYPE)*VOLUR(IEL)
         PV2=PERMY(ISTR,MTYPE)*VOLUR(IEL)
         PV3=PERMZ(ISTR,MTYPE)*VOLUR(IEL)
         PVK1=PV1*CKRWE(IEL)
         PVK2=PV2*CKRWE(IEL)
         PVK3=PV3*CKRWE(IEL)
         TPI=TETAF*PERMZ(ISTR,MTYPE)*IVOL(IEL)
         TVD=TETAF*VOLU(IEL)*RDT
         DO K=1,4
            KNOD=TETRA(K,IEL)
            PVB=PV1*BI(K,IEL)
            PVC=PV2*CI(K,IEL)
            PVD=PV3*DI(K,IEL)
            PVKB=PVK1*BI(K,IEL)
            PVKC=PVK2*CI(K,IEL)
            PVKD=PVK3*DI(K,IEL)
            SUM =0.0D0
            SUM1=0.0D0
            DO M=1,4
               MNOD=TETRA(M,IEL)
               SUM  = SUM + ( PVB*BI(M,IEL) + PVC*CI(M,IEL) + 
     1                        PVD*DI(M,IEL) ) * PTNEW(MNOD)
               SUM1 = SUM1 + LMASS(K,M)*(PNEW(MNOD) - PTIMEP(MNOD))
            END DO
            TSUMTD=TETAF*SUM + TPI*DI(K,IEL)
            SUM1TV=SUM1*TVD         
            DO L=1,4
               LNOD=TETRA(L,IEL)
               IND=TETJA(K,L,IEL)
               IF (LNOD .GE. KNOD) THEN
                  COEF1(IND)=COEF1(IND) + PVKB*BI(L,IEL) + 
     1                                    PVKC*CI(L,IEL) +
     2                                    PVKD*DI(L,IEL)
                  COEF2(IND)=COEF2(IND) + VE*LMASS(K,L)
               END IF
               COEF3(IND)=COEF3(IND) + DCKRW(LNOD)*TSUMTD + 
     1                                 DETAI(LNOD)*SUM1TV
            END DO
         END DO
         DO K=1,4
            KNOD=TETRA(K,IEL)
            DO L=1,4
               LNOD=TETRA(L,IEL)
               IF (LNOD .LT. KNOD) THEN
                  IND=TETJA(K,L,IEL)
                  INDS=TETJA(L,K,IEL)
                  COEF1(IND)=COEF1(INDS)
                  COEF2(IND)=COEF2(INDS)
               END IF
            END DO
         END DO
      END DO
C
      RETURN
      END
