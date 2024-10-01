C
C***********************************************************************
C             Subroutine ETRAN, computes root water uptake
C***********************************************************************
C
      SUBROUTINE ETRAN(N,NNOD,NSTR,ATMPOT,Z,PSI,PNODI,VEG_TYPE,QTRANIE)

      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INCLUDE 'SOILCHAR.H'

      INTEGER I,J,K,N,NNOD,NSTR
      INTEGER VEG_TYPE(*)
      REAL*8  SH2O,GX,ZERO,DZ,ZSURF,DEPTH,BETA
      REAL*8  S1,S2,GX1,GX2
      REAL*8  Z(*),PSI(*),PNODI(*)
      REAL*8  ATMPOT(*),QTRANIE(*)
      REAL*8  BTRAN(NNOD),BTRANI(N),ETP(NNOD),OMG(NNOD)
      DATA    ZERO/0.0d+00/
   
      CALL INIT0R(N,QTRANIE)
      CALL INIT0R(N,BTRANI)
      DO I = 1,NNOD
         BTRAN(I) = 0.0d0
         OMG(I)=0.0d0
         J = 1
         ZSURF = Z(I)
         DEPTH = 0.0d0
         IF (ATMPOT(I).LT.0.0d0) THEN
            ETP(I) = -1.0d0*SCF*ATMPOT(I)
         ELSE
            ETP(I) = 0.0d0
         END IF
         DO WHILE (DEPTH.LE.ZROOT(VEG_TYPE(I)))
            K = (J-1)*NNOD+I
            S1 = PCANA(VEG_TYPE(I))
            S2 = PCANA(VEG_TYPE(I))+1.0D-03
            IF (J.EQ.1) THEN
               DZ = (ZSURF-Z(K+NNOD))/2.0d0
            ELSE
               DZ = (Z(K-NNOD)-Z(K+NNOD))/2.0d0
            END IF
            SH2O = PSI(K)
            GX1 = (SH2O-PCWLT(VEG_TYPE(I))) / 
     &            (PCREF(VEG_TYPE(I))-PCWLT(VEG_TYPE(I)))
            GX1 = MIN(1.0d0,MAX(0.0d0,GX1))
            GX2 = 1.0d0 - (SH2O-S1)/(S2-S1)
            GX2 = MIN(1.0d0,MAX(0.0d0,GX2))
            GX = MIN(GX1,GX2)
CM          IF (SH2O.GE.PCANA) GX = 0.0d0
CM          PZ = -7.044D-03*ZROOT(I)**4 + 1.261D-01*ZROOT(I)**3
CM   +           -8.101D-01*ZROOT(I)**2 + 6.849D+00*ZROOT(I) - 3.345D+00
            BETA = (1-DEPTH/ZROOT(VEG_TYPE(I)))*
     &             DEXP(-1.0d0*PZ(VEG_TYPE(I))*DEPTH/ZROOT(VEG_TYPE(I)))
            BTRANI(K) = MAX(ZERO,BETA*DZ*GX)
            BTRAN(I)  = BTRAN(I) + BETA*DZ
            OMG(I)    = OMG(I) + GX*BETA*DZ
            J = J + 1
            IF ((J-1)*NNOD+I.LE.NMAX) THEN
               DEPTH = ZSURF - Z((J-1)*NNOD+I)
            ELSE
               write(6,*)'DECREASE ZROOT, TOO CLOSE TO THE BOTTOM'
               stop
            END IF
cm          write(111,*)k,gx1,gx2,beta,btrani(k)
         END DO
         BTRAN(I) = MAX(ZERO,BTRAN(I))
         OMG(I) = OMG(I)/BTRAN(I)
      END DO
      DO I=1,NNOD
cm       write(666,*)i,etp(i),btran(i),zroot(i)
         DO J=1,NSTR+1
            K = (J-1)*NNOD+I
            QTRANIE(K)=ETP(I)*BTRANI(K)/BTRAN(I)/
     &                 MAX(OMG(I),OMGC(VEG_TYPE(I)))
cm          write(111,*)k,qtranie(k),btrani(k),btran(i)
         END DO
      END DO

      RETURN
      END
