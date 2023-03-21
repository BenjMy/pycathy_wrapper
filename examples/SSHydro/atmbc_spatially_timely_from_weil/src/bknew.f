C
C**************************  BKNEW  ************************************
C
C  back-calculate fluxes at all Dirichlet nodes at the current time 
C  level. Newton case.
C
C***********************************************************************
C
      SUBROUTINE BKNEW(NP,NNOD,TETAF,TOPOL,JA,CONTP,IFATM,PDIFF,
     1                 QPNEW,QPOLD,ATMACT,ATMOLD,COEF1,XT5,
     2                 NSF,NSFNUM,NSFNOD,SFQ,SFQP,SFEX)
C
      IMPLICIT  NONE
      INCLUDE  'CATHY.H'
      INTEGER   I,J,K,M,II,JJ,INOD
      INTEGER   NP,NNOD,NSF
      INTEGER   TOPOL(*),JA(*),CONTP(*),IFATM(*)
      INTEGER   NSFNUM(*),NSFNOD(NSFMAX,*),SFEX(NSFMAX,*)
      REAL*8    TETAFR
      REAL*8    TETAF
      REAL*8    PDIFF(*),QPNEW(*),QPOLD(*),ATMACT(*),ATMOLD(*)
      REAL*8    COEF1(*),XT5(*)
      REAL*8    SFQ(NSFMAX,*),SFQP(NSFMAX,*)
C
      TETAFR=1.0D0/TETAF
      DO I=1,NP
         K=CONTP(I)
         QPNEW(I)=0.0D0
         II=TOPOL(K)
         JJ=TOPOL(K+1) - 1
         DO M=II,JJ
            QPNEW(I)=QPNEW(I) + COEF1(M)*PDIFF(JA(M))
         END DO
         QPNEW(I)=QPNEW(I) - XT5(K)
         QPNEW(I)=(QPNEW(I) - (1.0D0 - TETAF)*QPOLD(I))*TETAFR
      END DO
      DO I=1,NSF
         DO J=1,NSFNUM(I)
            IF (SFEX(I,J).EQ.1)THEN
               INOD=NSFNOD(I,J)
               SFQ(I,J)=0.0D0
               II=TOPOL(INOD)
               JJ=TOPOL(INOD+1) - 1
               DO M=II,JJ
                  SFQ(I,J)=SFQ(I,J) + COEF1(M)*PDIFF(JA(M))
               END DO
               SFQ(I,J)=SFQ(I,J) - XT5(INOD)
               SFQ(I,J)=(SFQ(I,J) - (1.0D0 - TETAF)*SFQP(I,J))*TETAFR
           END IF
         END DO
      END DO
      DO K=1,NNOD
         IF (IFATM(K) .EQ. 1 .OR. IFATM(K) .EQ. 2) THEN
            ATMACT(K)=0.0D0
            II=TOPOL(K)
            JJ=TOPOL(K+1) - 1
            DO M=II,JJ
               ATMACT(K)=ATMACT(K) + COEF1(M)*PDIFF(JA(M))
            END DO
            ATMACT(K)=ATMACT(K) - XT5(K)
            ATMACT(K)=(ATMACT(K) - (1.0D0 - TETAF)*ATMOLD(K))*TETAFR
         END IF
      END DO
C
      RETURN
      END
