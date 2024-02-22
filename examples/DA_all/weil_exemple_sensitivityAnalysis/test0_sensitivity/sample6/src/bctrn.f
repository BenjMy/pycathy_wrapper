C
C**************************  BCTRN  ************************************
C
C  impose Dirichlet and Neumann boundary conditions
C
C***********************************************************************
C
      SUBROUTINE BCTRN(nstep,NNOD,NPC,NNPC,NMC,NNMC,TOPOLC,
     1                 RMAX,TNOTIC,PC,MC1,COEF1C,JAC,ATMACT)
C
      IMPLICIT  NONE
      INTEGER   K,J,IND,nstep,i
      INTEGER   NPC,NMC,NNOD
      INTEGER   NNMC(*),NNPC(*),TOPOLC(*),JAC(*)
      REAL*8    RMAX
      REAL*8    TNOTIC(*),PC(*),MC1(*),COEF1C(*),ATMACT(*)
C
C  save a copy of RHS vector before imposing Dirichlet boundary
C  conditions (needed for back-calculation of fluxes)
C
C
C  Dirichlet BC's
C
      DO K=1,NPC
         J=NNPC(K)
         IND=TOPOLC(J)
         TNOTIC(J)=PC(K)*1.0D-9*RMAX
         COEF1C(IND)=1.0D-9*RMAX
c         write(189,*) 'k', k, ind
      END DO
C Neumann conditions
c     DO K=1,NMC
c           J=NNMC(K)

c           TNOTIC(J)=TNOTIC(J) + MC1(K)
c     END DO

C  Cauchy BC's (total flux in MC1, advective part in MC2).
C  Since the Cauchy values are taken to be zero at time 0, we
C  use an average flux value for the first time step.
C
      IF (NSTEP .EQ. 1) THEN
         DO K=1,NMC
            J=NNMC(K)
            TNOTIC(J)=TNOTIC(J) + 0.5D0*MC1(K)
         END DO
      ELSE
         DO K=1,NMC
            J=NNMC(K)
            TNOTIC(J)=TNOTIC(J) + MC1(K)
         END DO
      END IF
      
      DO K=1,NMC
         J=NNMC(K)
         IND=TOPOLC(J)-1
 300     IND=IND+1
         IF (JAC(IND) .NE. J) GO TO 300
         IF (IND .LE. NNOD) COEF1C(IND)= COEF1C(IND) - ATMACT(IND)
      END DO
C     
C     
      RETURN
      END
      
