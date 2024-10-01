C
C**************************  CONVER ************************************
C
C  check nonlinear convergence and switching of atmospheric
C  and seepage face boundary conditions
C
C***********************************************************************
C
      SUBROUTINE CONVER(ITER,N,NNOD,NSF,KSFCV,KSFCVT,NITER,ISFONE,
     1                  L2NORM,IKMAX,KSF,NSFNUM,NSFNOD,SFEX,SFEXIT,
     2                  SFFLAG,IFATM,SURF,PONDING,KSFZER,TOLSWI,
     3                  TIME,DELTAT,PIKMAX,PINF,PL2,FINF,FL2,PONDH_MIN,
     4                  PNEW,POLD,TNOTI,ATMACT,ATMPOT,SFQ,QTRANIE,
     5                  ARENOD,PONDNOD,OVFLNOD,DUPUIT,Z,PUNTDIRSF_NODE)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I
      INTEGER  ITER,N,NNOD,NSF,KSFCV,KSFCVT,NITER,ISFONE
      INTEGER  L2NORM,IKMAX,KSF,DUPUIT
      INTEGER  NSFNUM(*),NSFNOD(NSFMAX,*),PUNTDIRSF_NODE(*)
      INTEGER  SFEX(NSFMAX,*),SFEXIT(NSFMAX,*),SFFLAG(*),IFATM(*)
      LOGICAL  SURF,PONDING,KSFZER
      REAL*8   TOLSWI,TIME,DELTAT
      REAL*8   PIKMAX,PINF,PL2,FINF,FL2,PONDH_MIN
      REAL*8   PNEW(*),POLD(*),TNOTI(*),Z(*)
      REAL*8   ATMACT(*),ATMPOT(*),SFQ(NSFMAX,*),QTRANIE(*)
      REAL*8   ARENOD(*),PONDNOD(*),OVFLNOD(*)
      include 'NORMVL.H'
      include 'IOUNITS.H'
C    
C  calculate the nonlinear convergence and residual error norms
C     
      CALL NORMS(N,IKMAX,PIKMAX,PINF,PL2,FINF,FL2,PNEW,POLD,TNOTI)
      ITUMAX(ITER)=IKMAX
      PCURRV(ITER)=PNEW(IKMAX)
      PPREVV(ITER)=POLD(IKMAX)
      PIKMXV(ITER)=PIKMAX
      PL2V(ITER)=PL2
      FINFV(ITER)=FINF
      FL2V(ITER)=FL2
      WRITE(ITERM,1030) ITER,PL2,PIKMAX,IKMAX,PNEW(IKMAX),
     1                  POLD(IKMAX),FL2,FINF
      WRITE(IOUT4,1070) ITER,NITER,PL2,PINF,IKMAX,PNEW(IKMAX),
     1                  POLD(IKMAX),FL2,FINF
C
C  check for switching of atmospheric boundary conditions depending
C  on the setting of TOLSWI.
C  Note: SWITCH_OLD is the version of SWITCH for the uncoupled,
C  subsurface flow only version of the model. It should be superceded
C  by the new SWITCH, but for the time being we keep the old version
C  active as well.
C        
cd     write(99,*)'pnew prima di switch'
cd     do i=1,nnod
cd        write(99,*) i,pnew(i)
cd     end do
      IF ((L2NORM .EQ. 0 .AND. PINF .LE. TOLSWI)  .OR.
     1    (L2NORM .NE. 0 .AND. PL2  .LE. TOLSWI)) THEN 
cd        write(99,*) 'chiamata switch da conver'   
cd        write(99,*) 'time=',time,'iter=',iter
         IF (.NOT. SURF) THEN
cd            write(99,*) 'switch-old (not surf)'
            CALL SWITCH_OLD(NNOD,IFATM,ATMACT,ATMPOT,PNEW)
         ELSE
cd            write(99,*) 'switch-new (surf)' 
            CALL SWITCH(NNOD,IFATM,PONDING,TIME,DELTAT,PONDH_MIN,
     1                  ARENOD,PONDNOD,ATMPOT,ATMACT,QTRANIE,PNEW,
     2                  OVFLNOD)
         END IF
      END IF
C     
C  calculate new position of the exit point along each seepage
C  face and check for seepage face exit point convergence
C        
      IF (NSF .GT. 0) THEN
         IF (ISFONE .EQ. 0) THEN
            CALL EXTALL(N,NSF,NSFNUM,NSFNOD,SFEX,SFEXIT,SFQ,PNEW,
     1                  SFFLAG,ITER,TIME,DELTAT,KSF,DUPUIT,Z,
     2                  PUNTDIRSF_NODE)
         ELSE 
cm 
cm this subroutine is not used any longer, after the last update
cm of seepage face handling
cm
c           CALL EXTONE(N,NSF,NSFNUM,NSFNOD,SFEX,SFEXIT,SFQ,PNEW,
c    1                  SFFLAG,ITER,TIME,DELTAT,KSF,DUPUIT,Z,
c    2                  PUNTDIRSF_NODE)
            WRITE(*,*) 'CHANGE CODE IN conver.f FOR ISFONE=1'
            stop
         END IF
         IF (KSF .GT. 0) THEN  
            KSFCV=KSFCV + 1
            KSFCVT=KSFCVT + KSF
            KSFZER=.FALSE.
         ELSE     
            KSFZER=.TRUE.  
         END IF
      END IF
C
      RETURN
 1030 FORMAT(I6,2(1PE12.4),I6,2(1PE11.2),2(1PE11.3)) 
 1070 FORMAT(2I6,2(1PE10.3),I6,2(1PE10.2),2(1PE10.3))
      END
