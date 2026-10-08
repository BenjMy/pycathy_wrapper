C
C**************************  TRANSPMB **********************************
C
C calculate total in and out solute mass fluxes at the current time step
C
C***********************************************************************
C
      SUBROUTINE TRANSPMB(N,NNOD,NP,NQ,ACONTP,ACONTQ,QPNEW,Q,ATMACT,
     1                    IFATM,CNNEW,DELTAT,SFMASS,MASSIN,MASSOUT,
     2                    MASSINTOT,MASSOUTTOT)
C
      IMPLICIT  NONE
      INCLUDE  'CATHY.H'
      INTEGER   K
      INTEGER   N,NP,NQ,NNOD
      INTEGER   IFATM(*),ACONTP(*),ACONTQ(*)
      REAL*8    MASSIN,MASSOUT,MASSINTOT,MASSOUTTOT,DELTAT,SFMASS
      REAL*8    QPNEW(*),Q(*),CNNEW(*),ATMACT(*)
C
C  non-atmospheric Dirichlet nodes
C
      DO K=1,NP
         IF (QPNEW(K) .GT. 0.0D0) THEN
            MASSIN=MASSIN + QPNEW(K)*CNNEW(ACONTP(K))*DELTAT
         ELSE
            MASSOUT=MASSOUT + QPNEW(K)*CNNEW(ACONTP(K))*DELTAT
         END IF
      END DO
C
C  non-atmospheric Neumann nodes
C
      DO K=1,NQ
         IF (Q(K) .GT. 0.0D0) THEN
            MASSIN=MASSIN + Q(K)*CNNEW(ACONTQ(K))*DELTAT
         ELSE 
            MASSOUT=MASSOUT + Q(K)*CNNEW(ACONTQ(K))*DELTAT
         END IF
      END DO
C
C  atmospheric Dirichlet and Neumann nodes
C
      DO 100 K=1,NNOD
         IF (IFATM(K) .EQ. -1) GO TO 100
         IF (IFATM(K) .EQ. 1 .OR. IFATM(K) .EQ. 2) THEN
            IF (ATMACT(K) .GT. 0.0D0) THEN
               MASSIN=MASSIN + ATMACT(K)*CNNEW(K)*DELTAT
            ELSE
               MASSOUT=MASSOUT + ATMACT(K)*CNNEW(K)*DELTAT
            END IF
         ELSE
            IF (ATMACT(K) .GT. 0.0D0) THEN
               MASSIN=MASSIN + ATMACT(K)*CNNEW(K)*DELTAT
            ELSE
               MASSOUT=MASSOUT + ATMACT(K)*CNNEW(K)*DELTAT
            END IF
         END IF
  100 CONTINUE
C
      MASSOUT= MASSOUT + SFMASS
      MASSINTOT= MASSINTOT+ MASSIN
      MASSOUTTOT= MASSOUTTOT+ MASSOUT
C
      RETURN
      END
