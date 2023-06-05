C
C**************************  DETOUT ************************************
C
C  detailed output
C
C***********************************************************************
C
      SUBROUTINE DETOUT(KPRT,IPEAT,IPRT,N,NNOD,NUMVP,NSTR,NSTEP,TIME,
     1                 NODVP,SATSUR,PNEW,INDE,DEF,SW,CKRW,UNOD,VNOD,
     2                 WNOD,NT,UU,VV,WW,X,Y,Z,OVFLNOD,ATMACT,ARENOD,
     3                 PONDNOD,IFATM,PNODI,RECNOD,QTRANIE,CNNEW)
C
      IMPLICIT  NONE
      INCLUDE   'CATHY.H'
      INTEGER   I,K,KK,INOD
      INTEGER   KPRT,IPEAT,IPRT,N,NNOD,NUMVP,NSTR,NSTEP,NT
      INTEGER   NODVP(*),SATSUR(*),ifatm(*)
      REAL*8    TIME,FCDPORE
      REAL*8    PNEW(*),INDE(*),SW(*),CKRW(*),UNOD(*),VNOD(*),WNOD(*)
      REAL*8    DEF(*),PNODI(*),CNNEW(*)
      REAL*8    UU(*),VV(*),WW(*)
      REAL*8    X(*),Y(*),Z(*)
      real*8    ovflnod(*),atmact(*),arenod(*),pondnod(*),QTRANIE(*)
      real*8    RECNOD(NODMAX),ETA(NODMAX)
      INCLUDE  'IOUNITS.H'
C

      IF (IPRT .GE. 1) THEN
         WRITE(IOUT11,1000) NSTEP,TIME
         WRITE(IOUT11,1020) (PNEW(I),I=1,N)
         WRITE(IOUT42,1000) NSTEP,TIME
         WRITE(IOUT42,1020) (PONDNOD(I),I=1,NNOD)
         IF (IPRT .GE. 2) THEN
            WRITE(IOUT12,1000) NSTEP,TIME
            WRITE(IOUT12,1040) (UNOD(I),VNOD(I),WNOD(I),I=1,N)
            IF (IPRT .GE. 3) THEN 
ccc               WRITE(IOUT14,1000) NSTEP,TIME
ccc               WRITE(IOUT14,1020) (CKRW(I),I=1,N)
               IF (IPRT .GE. 4) THEN
                  WRITE(IOUT15,1010) TIME 
                  WRITE(IOUT15,1020) (UU(I),I=1,NT)
                  WRITE(IOUT15,1020) (VV(I),I=1,NT)
                  WRITE(IOUT15,1020) (WW(I),I=1,NT)
                  WRITE(IOUT13,1000) NSTEP,TIME
                  WRITE(IOUT13,1020) (SW(I),I=1,N)
               END IF
            END IF
         END IF
      END IF
      IF (IPEAT .EQ. 1)  THEN
         WRITE(IOUTPT,1000) NSTEP,TIME
         WRITE(IOUTPT,1020) (INDE(I),I=1,N)
      END IF
      IF (NUMVP .GT. 0 .AND. IPEAT .EQ. 0) THEN
         WRITE(IOUT6,1000) NSTEP,TIME
         DO I=1,NUMVP
            INOD=NODVP(I)
            IF (INOD .GE. 1 .AND. INOD .LE. NNOD) THEN
               WRITE(IOUT6,1070) INOD,X(INOD),Y(INOD)
               DO K=0,NSTR
                  KK=K*NNOD + INOD
                  WRITE(IOUT6,1080) Z(KK),PNEW(KK),SW(KK),CKRW(KK),
     +                              QTRANIE(KK),CNNEW(KK)
C                 WRITE(100+KPRT,1080) Z(KK),PNEW(KK),SW(KK),CKRW(KK)
               END DO
            END IF
         END DO
      END IF
      IF (NUMVP .GT. 0 .AND. IPEAT .EQ. 1) THEN
         WRITE(IOUT6,1000) NSTEP,TIME
         DO I=1,NUMVP
            INOD=NODVP(I)
            IF (INOD .GE. 1 .AND. INOD .LE. NNOD) THEN
               WRITE(IOUT6,1075) INOD,X(INOD),Y(INOD)
               DO K=0,NSTR
                  KK=K*NNOD + INOD
                  WRITE(IOUT6,1090) Z(KK),PNEW(KK),SW(KK),CKRW(KK),
     1                              INDE(KK),DEF(KK),(SW(KK)*
     2                  FCDPORE(INDE(KK),SW(KK),(Z(INOD)-Z(KK)),
     3                  PNODI(KK)))
               END DO
            END IF
         END DO
      END IF
ccc   write(55,*) 'time=',time
ccc   write(66,*) 'time=',time
ccc   write(77,*) 'time=',time
ccc   write(88,*) 'time=',time
ccc   do i=1,nnod
ccc      write(55,*)i,ovflnod(i)/arenod(i)
ccc      write(66,*)i,atmact(i)/arenod(i)
ccc      write(77,*) i,pondnod(i)
ccc      write(88,*) i,ifatm(i)
ccc   end do
      WRITE(IOUT16,1000) NSTEP,TIME
      WRITE(IOUT17,1000) NSTEP,TIME
      WRITE(IOUT18,1000) NSTEP,TIME
      WRITE(IOUT44,1000) NSTEP,TIME
      WRITE(777,1000) NSTEP,TIME
      WRITE(IOUT16,2000) 
      WRITE(IOUT17,2020) 
      WRITE(IOUT18,2040)
      WRITE(IOUT44,2041)
      WRITE(777,2042)
      DO I=1,NNOD
         SATSUR(I)=1
         IF (PNEW(I) .GE. 0.0D0) THEN
            SATSUR(I)=3
            DO K=1,NSTR
               KK=K*NNOD + I
               IF (PNEW(KK) .LT. 0.0D0) SATSUR(I)=2
            END DO
         END IF
         WRITE(IOUT16,2060) I,X(I),Y(I),PNEW(I)
         WRITE(IOUT17,2080) I,X(I),Y(I),SATSUR(I)
         WRITE(IOUT18,2060) I,X(I),Y(I),SW(I)
         WRITE(IOUT44,2060) I,X(I),Y(I),RECNOD(I)/ARENOD(I)
         ETA(I)=0.0D0
         DO K=1,NSTR+1
            ETA(I)=ETA(I)+QTRANIE((K-1)*NNOD+I)
         END DO
         WRITE(777,2060) I,X(I),Y(I),ETA(I)/ARENOD(I)
      END DO
C
      RETURN
 1000 FORMAT(I7,1PE16.8,'     NSTEP   TIME')
 1010 FORMAT(1PE16.8,'       TIME')
 1020 FORMAT(5(1PE15.6))
 1040 FORMAT(3(1PE15.6))
 1070 FORMAT(' SURFACE NODE = ',I5,'  X = ',1PE12.4,'  Y = ',1PE12.4,/,
     1       '              Z  PRESSURE HEAD             SW',
     2       '           CKRW        QTRANIE')
 1075 FORMAT(' SURFACE NODE = ',I5,'  X = ',1PE12.4,'  Y = ',1PE12.4,/,
     1       '              Z  PRESSURE HEAD             SW',
     2       '           CKRW           INDE            DEF',
     3       '    STOR1/STOR2')
 1080 FORMAT(6(1PE15.6))
 1090 FORMAT(7(1PE15.6))
 2000 FORMAT(' SURFACE NODE              X              Y',
     1       '  PRESSURE HEAD')
 2020 FORMAT(' SURFACE NODE              X              Y',
     1       '         SATSUR')
 2040 FORMAT(' SURFACE NODE              X              Y',
     1       '             SW')
 2041 FORMAT(' SURFACE NODE              X              Y',
     1       '      REC. FLUX')
 2042 FORMAT(' SURFACE NODE              X              Y',
     1       '      ACT. ETRA')
 2050 FORMAT(' SURFACE NODE              X              Y',
     1       '  PONDING HEAD ')
 2060 FORMAT(7X,I6,3(1PE15.6),i6,1pe15.6)
 2080 FORMAT(7X,I6,2(1PE15.6),13X,I2)
      END
