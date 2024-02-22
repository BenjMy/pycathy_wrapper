C
C**************************  FLUXMB ************************************
C
C  calculate total inflow and outflow fluxes at the current time level
C
C***********************************************************************
C
      SUBROUTINE FLUXMB(N,NP,NQ,NNOD,NSF,NUDC,QPNEW,Q,IFATM,ATMACT,
     1                  NDIN,NDOUT,NNIN,NNOUT,ADIN,ADOUT,ANIN,ANOUT,
     2                  NUDIN,NUDOUT,NUDNOD,NUDCUM,
     3                  NSFNUM,NSFNOD,SFEX,SFQ,SFFLW,SFFLAG,
     4                  TIME,DELTAT)
C
      IMPLICIT  NONE
      INCLUDE  'CATHY.H'
      INTEGER   I,J,K
      INTEGER   N,NP,NQ,NNOD,NSF,NUDC
      INTEGER   IFATM(*),NSFNUM(*),NSFNOD(NSFMAX,*),SFEX(NSFMAX,*)
      INTEGER   SFFLAG(*)
      REAL*8    NDIN,NDOUT,NNIN,NNOUT,ADIN,ADOUT,ANIN,ANOUT
      REAL*8    NUDIN,NUDOUT
      REAL*8    SFFLW,TIME,DELTAT
      REAL*8    QPNEW(*),Q(*),ATMACT(*)
      REAL*8    NUDNOD(*),NUDCUM(*),SFQ(NSFMAX,*)
      INCLUDE  'IOUNITS.H'
C
C  non-atmospheric Dirichlet nodes
C
      NDIN=0.0D0
      NDOUT=0.0D0
      DO K=1,NP
         IF (QPNEW(K) .GT. 0.0D0) THEN
            NDIN=NDIN + QPNEW(K)
         ELSE
            NDOUT=NDOUT + QPNEW(K)
         END IF
      END DO
C
C  non-atmospheric Neumann nodes
C
      NNIN=0.0D0
      NNOUT=0.0D0
      DO K=1,NQ
         IF (Q(K) .GT. 0.0D0) THEN
            NNIN=NNIN + Q(K)
         ELSE 
            NNOUT=NNOUT + Q(K)
         END IF
      END DO
C
C  seepage face nodes (only actual seepage face nodes; potential
C  seepage face nodes have zero flux)
C
      SFFLW=0.0D0
      DO I=1,NSF
         DO J=1,NSFNUM(I)
            IF(SFEX(I,J).EQ.1) THEN 
               SFFLW=SFFLW + SFQ(I,J)
               IF (SFQ(I,J) .GE. 0.0D0) THEN
               SFFLAG(5)=SFFLAG(5) + 1
c               WRITE(IOUT10,2500) I,TIME,DELTAT,SFQ(I,J),J,NSFNOD(I,J)
               END IF
            END IF
         END DO
      END DO
C
C  atmospheric Dirichlet and Neumann nodes
C
      ADIN=0.0D0
      ADOUT=0.0D0
      ANIN=0.0D0
      ANOUT=0.0D0
      DO 100 K=1,NNOD
         IF (IFATM(K) .EQ. -1) GO TO 100
         IF (IFATM(K) .EQ. 1 .OR. IFATM(K) .EQ. 2) THEN
            IF (ATMACT(K) .GT. 0.0D0) THEN
               ADIN=ADIN + ATMACT(K)
            ELSE 
               ADOUT=ADOUT + ATMACT(K)
            END IF
         ELSE
            IF (ATMACT(K) .GT. 0.0D0) THEN
               ANIN=ANIN + ATMACT(K)
            ELSE
               ANOUT=ANOUT + ATMACT(K)
            END IF
         END IF
  100 CONTINUE
C
C  nudging term contribution
C
      NUDIN=0.0D0
      NUDOUT=0.0D0
      IF (NUDC .GT. 0) THEN
         DO K=1,N
            IF (NUDNOD(K) .GT. 0.0D0) THEN
               NUDIN=NUDIN + NUDNOD(K)
            ELSE 
               NUDOUT=NUDOUT + NUDNOD(K)
            END IF
            NUDCUM(K)=NUDCUM(K) + NUDNOD(K)
         END DO
      END IF
C
      RETURN
 2500 FORMAT(  ' SFFLAG(5) AT SEEPAGE FACE ',I4,
     1         ' (TIME=',1PE8.2,', DELTAT=',1PE8.2,')',
     2       /,11X,'NON-NEGATIVE FLUX (INFLOW) OF ',1PE9.3,' AT NODE ',
     3         I3,' (NODE # ',I6,')')
      END
