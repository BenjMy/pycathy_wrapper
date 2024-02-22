C
C**************************  BCPIC  ************************************
C
C  impose boundary conditions, and
C  save diagonal elements of LHS system matrix corresponding to
C  Dirichlet nodes: Picard scheme
C
C***********************************************************************
C
      SUBROUTINE BCPIC(NP,NQ,CONTP,CONTQ,LHSP,Q,NNOD,N,
     1                 LHSATM,IFATM,ATMACT,ATMOLD,TOPOL,COEF1,TNOTI,
     2                 TETAF,RMAX,NUMDIR,NODDIR,
     3                 NSF,NSFNUM,NSFNOD,SFEX,LHSSF,QTRANIE)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I,J,K,IND,INOD
      INTEGER  NP,NQ,NNOD,NUMDIR,NSF,N
      INTEGER  CONTP(*),CONTQ(*),IFATM(*),TOPOL(*),NODDIR(*)
      INTEGER  NSFNUM(*),NSFNOD(NSFMAX,*),SFEX(NSFMAX,*)
      REAL*8   TETAF,RMAX
      REAL*8   LHSP(*),Q(*),LHSATM(*),ATMACT(*),ATMOLD(*)
      REAL*8   COEF1(*),TNOTI(*),LHSSF(NSFMAX,*),QTRANIE(NMAX)
C
C  condizioni di Dirichlet: assegna al nodo j-esimo il valore 
C  nullo della pressione e lo sostitisce al termine noto.
C  moltiplicando il coeff. diagonale della j-esima riga
C  per un numero molto grande (1.0e-9*RMAX) si trasforma
C  praticamente l'equazione j-esima nell'identita': PDIFF(j)=0
C  Note: since we solve for pressure head differences,
C  the Dirichlet BC's imposed have a value of zero.
C
      DO K=1,NP
         J=CONTP(K)
         IND=TOPOL(J)
         TNOTI(J)=0.0D0
         LHSP(K)=COEF1(IND)
         COEF1(IND)=1.0D-9*RMAX
      END DO
      NUMDIR=NP
      DO I=1,NSF
         DO J=1,NSFNUM(I)
            IF (SFEX(I,J).EQ.1) THEN
               INOD=NSFNOD(I,J)
               NUMDIR=NUMDIR + 1
               NODDIR(NUMDIR)=INOD
               IND=TOPOL(INOD)
               TNOTI(INOD)=0.0D0
               LHSSF(I,J)=COEF1(IND)
               COEF1(IND)=1.0D-9*RMAX
            END IF
         END DO
      END DO
      DO K=1,NNOD
         IF (IFATM(K) .EQ. 1 .OR. IFATM(K) .EQ. 2) THEN
            NUMDIR=NUMDIR + 1
            NODDIR(NUMDIR)=K
            IND=TOPOL(K)
            TNOTI(K)=0.0D0
            LHSATM(K)=COEF1(IND)
            COEF1(IND)=1.0D-9*RMAX
         END IF
      END DO
C
C  condizioni di Neumann: si assegna al nodo j la portata ottenuta 
C  come media tra le portate assegnate al tempo t e t+dt e si somma 
C  al corrispondente termine noto
C  Note: no need to do this for potential seepage face nodes since
C  the assigned Neumann flux is 0.0.
C
      DO K=1,NQ
         J=CONTQ(K)
         TNOTI(J)=TNOTI(J) + Q(K)
      END DO
      DO K=1,NNOD
         IF (IFATM(K) .EQ. 0) THEN
               TNOTI(K)=TNOTI(K) + 
     1                  (TETAF*ATMACT(K) + (1.0D0 - TETAF)*ATMOLD(K))
         END IF
      END DO
C
C  root extraction of water through leaf stomata
C
      DO K=1,N
          TNOTI(K) = TNOTI(K) - QTRANIE(K)
      END DO
C
      RETURN
      END
