C
C**************************  NUDPIC ************************************
C
C  add nudging contribution to the RHS vector: Picard scheme
C
C***********************************************************************
C
      SUBROUTINE NUDPIC(N,NT,NUDN,NUDC,NUDCTR,TIME,NUDG,NUDFLAG,WFLAG,
     1                  TETRA,NUDTET,IVOL,PTNEW,
     2                  X,Y,Z,VOLUR,VOLU,SW,PNODI,TNOTI,AI,BI,CI,DI,
     3                  NUDTIM,NUDRXY,NUDRZ,NUDX,NUDY,NUDZ,NUDEPS,
     4                  NUDDIF,NUDSMC,NUDNOD,NUDVAL,NUDTAU)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I,IN,INOD,K,KE
      INTEGER  N,NT,NUDN,NUDC,NUDCTR,NUDFLAG,WFLAG
      INTEGER  TETRA(5,*),NUDTET(*),IVOL(*)
      REAL*8   NUDWXY,NUDWZ,NUDWT
      REAL*8   NUDNORM,NUDELT,NUDGE_WINK
      REAL*8   TIME,NUDG
      REAL*8   PTNEW(*),X(*),Y(*),Z(*),VOLUR(*),VOLU(*)
      REAL*8   SW(*),PNODI(*),TNOTI(*)
      REAL*8   AI(4,*),BI(4,*),CI(4,*),DI(4,*)
      REAL*8   NUDTIM(*),NUDRXY(*),NUDRZ(*)
      REAL*8   NUDX(*),NUDY(*),NUDZ(*),NUDEPS(*)
      REAL*8   NUDDIF(*),NUDSMC(*),NUDNOD(*)
      REAL*8   NUDVAL(MAXNUDC,*),NUDTAU(MAXNUDT,*)
C
      IF (NUDC .EQ. 0) GO TO 900
      CALL INIT0R(N,NUDNOD)
C
C  evaluate model (computed) results at the observation points
C
      CALL NUDCPT(NUDFLAG,NUDN,TETRA,NUDTET,IVOL,VOLUR,SW,PNODI,
     1            PTNEW,AI,BI,CI,DI,NUDX,NUDY,NUDZ,NUDSMC)
C
C  compute the nudging term
C
      DO K=1,NUDC
         DO I=1,NUDN
CM          IF (NUDFLAG. EQ. 0) THEN
               NUDDIF(I)=NUDG * NUDEPS(I) * (NUDVAL(K,I) - NUDSMC(I))
CM          ELSE
CM             NUDDIF(I)=NUDG * NUDEPS(I) * (NUDVAL(K,I) - NUDSMC(I)) *
CM   1                   NUDETA(I)
CM          END IF
         END DO
         DO KE=1,NT
            NUDNORM=0.0D0
            NUDELT=0.0D0
            DO I=1,NUDN
               DO IN=1,4
                  INOD=TETRA(IN,KE)
                  NUDGE_WINK=NUDWXY(WFLAG,X(INOD),Y(INOD),NUDX(I),
     1                              NUDY(I),NUDRXY(I)) *
     2                       NUDWZ(WFLAG,Z(INOD),NUDZ(I),NUDRZ(I)) *
     3                       NUDWT(WFLAG,TIME,NUDTIM(NUDCTR-NUDC+K),
     4                             NUDTAU(NUDCTR-NUDC+K,I))
                  NUDNORM=NUDNORM + NUDGE_WINK
                  NUDELT=NUDELT + NUDGE_WINK*NUDGE_WINK*NUDDIF(I)
               END DO
            END DO
            IF (NUDNORM .NE. 0.0D0) THEN
               NUDELT=NUDELT / NUDNORM
            ELSE 
               NUDELT=0.0D0
            END IF
            DO IN=1,4
               INOD=TETRA(IN,KE)
               NUDNOD(INOD)=NUDNOD(INOD) + NUDELT*VOLU(KE)*0.25D0
            END DO
         END DO
      END DO
C
C  add the nudging term to the RHS vector
C
      DO I=1,N
         TNOTI(I)=TNOTI(I) + NUDNOD(I)
      END DO
C
  900 RETURN
      END
