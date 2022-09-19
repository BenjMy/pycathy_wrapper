C
C**************************  SHLNEW ************************************
C
C  restore diagonal elements of LHS system matrix corresponding to
C  Dirichlet nodes. 
C  Also, set solution at Dirichlet nodes to the prescribed 
C  values. This is done since the solved solution at
C  Dirichlet nodes may not be exactly equal to the presribed
C  values, due to inaccuracies and roundoff errors which could 
C  arise from the way we imposed Dirichlet conditions (by multiplying
C  diagonal terms by a 'large' number).
C  Newton scheme (nonsymmetric storage).
C
C***********************************************************************
C
      SUBROUTINE SHLNEW(NP,TOPOL,JA,CONTP,PRESC,LHSP,COEF1,PNEW,POLD,
     1                  NSF,NSFNUM,NSFNOD,SFEX,LHSSF,DUPUIT,Z,
     2                  NNOD,IFATM,LHSATM)
C
      IMPLICIT  NONE
      INCLUDE  'CATHY.H'
      INTEGER   I,J,IND,INOD
      INTEGER   NP,NSF,NNOD,DUPUIT
      INTEGER   TOPOL(*),JA(*),CONTP(*),NSFNUM(*),NSFNOD(NSFMAX,*)
      INTEGER   SFEX(NSFMAX,*),IFATM(*)
      REAL*8    PRESC(*),LHSP(*),COEF1(*),PNEW(*),POLD(*),Z(*)
      REAL*8    LHSSF(NSFMAX,*),LHSATM(*)
C
      DO I=1,NP
         J=CONTP(I)
         PNEW(J)=PRESC(I)
         IND=TOPOL(J) - 1
  300    IND=IND + 1
         IF (JA(IND) .NE. J) GO TO 300
         COEF1(IND)=LHSP(I)
      END DO
      DO I=1,NSF
         DO J=1,NSFNUM(I)
            IF (SFEX(I,J).EQ.1)THEN
               INOD=NSFNOD(I,J)
               IF (DUPUIT.EQ.0) THEN
                  PNEW(INOD)=0.0D0
c              ELSE
c                PNEW(INOD)=0.0D0+Z(NSFNOD(I,SFEX(I)))-Z(INOD)
               END IF
               IND=TOPOL(INOD) - 1
  350          IND=IND + 1
               IF (JA(IND) .NE. INOD) GO TO 350
               COEF1(IND)=LHSSF(I,J)
            END IF
         END DO
      END DO
      DO I=1,NNOD
         IF (IFATM(I) .EQ. 1 .OR. IFATM(I) .EQ. 2) THEN
            PNEW(I)=POLD(I)
            IND=TOPOL(I) - 1
  400       IND=IND + 1
            IF (JA(IND) .NE. I) GO TO 400
            COEF1(IND)=LHSATM(I)
         END IF
      END DO
C
      RETURN
      END
