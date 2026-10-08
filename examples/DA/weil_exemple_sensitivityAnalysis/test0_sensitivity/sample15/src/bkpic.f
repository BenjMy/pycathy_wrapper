C
C**************************  BKPIC  ************************************
C
C  back-calculate fluxes at all Dirichlet nodes at the current time 
C  level. Picard case.
C
C***********************************************************************
C
      SUBROUTINE BKPIC(N,NP,NNOD,TETAF,TOPOL,JA,CONTP,IFATM,PDIFF,
     1                 QPNEW,QPOLD,ATMACT,ATMOLD,COEF1,XT5,SCR,
     2                 NSF,NSFNUM,NSFNOD,SFQ,SFQP,SFEX)
C
      IMPLICIT  NONE
      INCLUDE  'CATHY.H'
      INTEGER   I,J,K,M,JM,INOD
      INTEGER   N,NP,NNOD,NSF
      INTEGER   TOPOL(*),JA(*),CONTP(*),IFATM(*)
      INTEGER   NSFNUM(*),NSFNOD(NSFMAX,*),SFEX(NSFMAX,*)
      REAL*8    TETAFR
      REAL*8    TETAF
      REAL*8    PDIFF(*),QPNEW(*),QPOLD(*),ATMACT(*),ATMOLD(*)
      REAL*8    COEF1(*),XT5(*),SCR(*)
      REAL*8    SFQ(NSFMAX,*),SFQP(NSFMAX,*)
C
      TETAFR=1.0D0/TETAF
      CALL INIT0R(N,SCR)
      DO K=1,N
         I=TOPOL(K)
         J=TOPOL(K+1) - 1
         DO M=I,J
            JM=JA(M)
            SCR(K)=SCR(K) + COEF1(M)*PDIFF(JM)
            IF (M .NE. I) SCR(JM)=SCR(JM) + COEF1(M)*PDIFF(K)
         END DO
         SCR(K)=SCR(K) - XT5(K)
      END DO
      DO I=1,NP
         K=CONTP(I)
         QPNEW(I)=(SCR(K) - (1.0D0 - TETAF)*QPOLD(I))*TETAFR
      END DO
      DO I=1,NSF
         DO J=1,NSFNUM(I)
            IF (SFEX(I,J).EQ.1) THEN
               INOD=NSFNOD(I,J)
               SFQ(I,J)=(SCR(INOD) - (1.0D0 - TETAF)*SFQP(I,J))*TETAFR
            END IF
         END DO
      END DO
      DO K=1,NNOD
         IF (IFATM(K) .EQ. 1 .OR. IFATM(K) .EQ. 2) THEN
             ATMACT(K)=(SCR(K) - (1.0D0 - TETAF)*ATMOLD(K))*TETAFR
         END IF
      END DO
C
      RETURN
      END
