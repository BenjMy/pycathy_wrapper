        subroutine ut3d(iout6,NNOD,nstep,n,ntri3,nt,ntermc,
     1             npc,nmc,
     1             niaux,nraux,tetra,ip4,iaux,iparm,aux,
     2             tetjac,jac,nnpc,nnmc,topolc,
     3             diffus,gammas,rmax,deltat,tetac,tolcg,
     4             x,y,z,uu,vv,ww,swe,volu,volur,volnod,ccur,
     5             ctimep,time,xt5c,tnotic,pc,mc1,poros,
     6             kd,alfal,alfat,lmassc,lhsc,coef1c,
     7             coef2c,cpuvec,iprt1,supdiffus,atmact)

C
C    Finite Element for the dispersive part of the transport equation
C
      implicit none
      include 'CATHY.H'
      INTEGER   iscratch(1),I,niaux,nraux,iout6,j
      integer   iparm(10), IAUX(*),ex_num
      INTEGER   N,NT,NTERMC,NPC,NMC,NNOD
      INTEGER   IPRT1,nstep,ntri3
      INTEGER   TETRA(5,*),IP4(4,4),TETJAC(4,4,*),JAC(*)
      INTEGER   NNPC(*),NNMC(*),TOPOLC(*)
      REAL      CPUVEC(*),CPUT1
      REAL*8    DIFFUS,GAMMAS,RMAX,DELTAT,TETAC,TOLCG
      REAL*8    X(*),Y(*),Z(*),UU(*),VV(*),WW(*)
      real*8    SWE(*),VOLU(*),VOLUR(*),volnod(*)
      REAL*8    CCUR(*),CTIMEP(*),XT5C(*),TNOTIC(*),ATMACT(*)
cxcx  REAL*8    PC(*),MC1(*),aux(*),supdiffus,prec(MAXTRM)
      REAL*8    PC(*),MC1(*),aux(*),supdiffus
      REAL*8    POROS(MAXSTR,*),KD(MAXSTR,*)
      REAL*8    ALFAL(MAXSTR,*),ALFAT(MAXSTR,*)
      REAL*8    LMASSC(4,4),LHSC(*),COEF1C(*),COEF2C(*)
      real*8    time, zero
      
      zero=0.d0

C  initialize the arrays which will contain the global stiffness
C  and mass matrices
C
      CALL TIM(CPUT1,1)
      CALL INIT0(NTERMC,COEF1C)
      CALL INIT0(NTERMC,COEF2C)
      CALL TIM(CPUT1,2)
      CPUVEC(1)=CPUVEC(1) + CPUT1
C
C  assemble global stiffness and mass matrices from the local
C  contributions of each element
C
      CALL TIM(CPUT1,1)
      CALL ASSUT3(NT,NTRI3,TETRA,IP4,TETJAC,DIFFUS,GAMMAS,
     1            X,Y,Z,UU,VV,WW,SWE,VOLU,VOLUR,COEF1C,
     2            COEF2C,POROS,KD,ALFAL,ALFAT,LMASSC,supdiffus)
      CALL TIM(CPUT1,2)
      CPUVEC(2)=CPUVEC(2) + CPUT1

C  assemble RHS vector, without the contributions of the boundary
C  conditions
C
      CALL TIM(CPUT1,1)
      CALL RHSTRN(N,TOPOLC,JAC,DELTAT,TETAC,TNOTIC,CTIMEP,
     1            COEF1C,COEF2C)
      CALL TIM(CPUT1,2)
      CPUVEC(3)=CPUVEC(3) + CPUT1
C
C  assemble the global LHS system matrix from the stiffness and mass
C  matrices
      CALL TIM(CPUT1,1)
      CALL CFMAT(NTERMC,TETAC,DELTAT,COEF1C,COEF2C)
      CALL TIM(CPUT1,2)
      CPUVEC(4)=CPUVEC(4) + CPUT1
C
C  save diagonal elements of LHS system matrix corresponding to
C  Dirichlet nodes
C
      CALL TIM(CPUT1,1)
      CALL LHSTRN(NPC,NNPC,TOPOLC,JAC,LHSC,COEF1C)
C
C  impose boundary conditions
C
      CALL BCTRN(nstep,NNOD,NPC,NNPC,NMC,NNMC,TOPOLC,
     1           RMAX,TNOTIC,PC,MC1,COEF1C,JAC,ATMACT)
      CALL TIM(CPUT1,2)
      CPUVEC(5)=CPUVEC(5) + CPUT1
      CALL TIM(CPUT1,1)
c----------------------

c
c * PCG solution of the linear system
c * we have to pass iscratch(1)=0 and the vector iscratch to SYMSLV
c * to ensure that no elements are excluded in the calculation of the
c * norm of the residual
c
         iscratch(1)=0
c     con cgsolv i risultati sono errati... perche'? quali parametri
c     errati?s
      call cgsolv(n,ntermc,0,niaux,nraux,iparm,iscratch,
     1     topolc,jac,iaux,tolcg,aux,coef1c,tnotic,ccur)
c      call kersh(100,n,ntermc,topolc,jac,coef1c,prec)
c      call grad(maxnod,n,ntermc,2000,topolc,jac,coef1c,prec,tnotic,
c     1       tolcg,ccur,rmax)
c       do j=1,n
c          write(190,*) ccur(j)
c       end do
C
C  restore diagonal elements of LHS system matrix corresponding to
C  Dirichlet nodes.
C  Also, set solution at Dirichlet nodes to the prescribed
C  values. This is done since the solved solution at
C  Dirichlet nodes may not be exactly equal to the presribed
C  values, due to inaccuracies and roundoff errors which could
C  arise from the way we imposed Dirichlet conditions (by multiplying
C  diagonal and RHS terms by a 'large' number).
C
      CALL SHLSYM(NPC,TOPOLC,NNPC,PC,LHSC,COEF1C,CCUR)
      CALL TIM(CPUT1,2)
      CPUVEC(6)=CPUVEC(6) + CPUT1
C
      RETURN
      END

