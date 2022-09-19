       subroutine ASSUT3(NT,NTRI3,TETRA,IP4,TETJAC,DIFFUS,GAMMAS,
     1            X,Y,Z,UU,VV,WW,SWE,VOLU,VOLUR,COEF1C,
     2            COEF2C,POROS,KD,ALFAL,ALFAT,LMASSC,supdiffus)
       IMPLICIT NONE
       include 'CATHY.H'
       INTEGER   IEL,I,K,L,MTYPE,KNOD,LNOD,IND,INDS
       INTEGER NT, ISTR, NTRI3,IR
       INTEGER   TETRAL(5)
       INTEGER   TETRA(5,*),IP4(4,4),TETJAC(4,4,*)
       REAL*8    APR,APRL,DXXA,DYYA,DZZA,ASIGN3,DXXAB,DYYAB,DZZAB
       REAL*8    VEL,DIFFUS,GAMMAS,supdiffus
       REAL*8    XX(4),YY(4),ZZ(4),BIL(4),CIL(4),DIL(4)
       REAL*8    X(*),Y(*),Z(*),UU(*),VV(*),WW(*),SWE(*)
       REAL*8    VOLU(*),VOLUR(*)
       REAL*8    COEF1C(*),COEF2C(*)
       REAL*8    POROS(MAXSTR,*),KD(MAXSTR,*)
       REAL*8    ALFAL(MAXSTR,*),ALFAT(MAXSTR,*)
       REAL*8    LMASSC(4,4)

       DO IEL=1,NT
         MTYPE=TETRA(5,IEL)
         ISTR=1 + IEL/NTRI3
         IR=MOD(IEL,NTRI3)
         IF (IR .EQ. 0) ISTR=ISTR-1
         DO I=1,4
            XX(I)=X(TETRA(I,IEL))
            YY(I)=Y(TETRA(I,IEL))
            ZZ(I)=Z(TETRA(I,IEL))
            TETRAL(I)=I
         END DO
C  rotate the current element so that the local x-axis is aligned with
C  the velocity vector, and re-compute the basis function
C  coefficients (divided by 6) for this rotated element
C
         CALL ANIS3D(UU(IEL),VV(IEL),WW(IEL),XX,YY,ZZ,VEL)
         CALL BASIS3D(IP4,TETRAL,XX,YY,ZZ,BIL,CIL,DIL)
C
C  assemble the mass matrix and the stiffness matrix
C
         APR=VOLU(IEL)*(POROS(ISTR,MTYPE)*SWE(IEL) + GAMMAS*
     1       (1.0D0 - POROS(ISTR,MTYPE))*KD(ISTR,MTYPE))
         DXXA=(VEL*ALFAL(ISTR,MTYPE) + DIFFUS*POROS(ISTR,MTYPE)*
     1        SWE(IEL))*VOLUR(IEL)
         DYYA=(VEL*ALFAT(ISTR,MTYPE) + DIFFUS*POROS(ISTR,MTYPE)*
     1        SWE(IEL))*VOLUR(IEL)
         DZZA=DYYA
         supdiffus= max(max(dabs(dxxa),dabs(dzza)),supdiffus)
         DO K=1,4
            KNOD=TETRA(K,IEL)
            DXXAB= DXXA*BIL(K)
            DYYAB= DYYA*CIL(K)
            DZZAB= DZZA*DIL(K)
            DO L=1,4
               LNOD=TETRA(L,IEL)
               IF (LNOD .GE. KNOD) THEN
               IND=TETJAC(K,L,IEL)
               COEF1C(IND)=COEF1C(IND)+DXXAB*BIL(L)+DYYAB*CIL(L)
     1                     + DZZAB*DIL(L)
               COEF2C(IND)=COEF2C(IND)+ APR*LMASSC(K,L)
               
               END IF
            END DO
         END DO
       END DO
c       write(*,*) coef2C(300)
c       write(*,*) supdiffus
       RETURN
       END



