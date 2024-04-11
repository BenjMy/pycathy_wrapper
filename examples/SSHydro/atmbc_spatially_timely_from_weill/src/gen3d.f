C
C**************************  GEN3D *************************************
C
C  la routine GEN3D crea il file di topologia per gli elementi 
C  tetraedrici di un acquifero partendo dalla topologia dei triangoli 
C  di base registrata col seguente ordine
C     1 punto: punto di convergenza della coppia di
C              frecce aventi rotazione opposta
C     2 punto: punto di convergenza della 3m freccia
C
C***********************************************************************
C
      SUBROUTINE GEN3D(N,NT,NNOD,NSTR,NTRI,IPRT1,IVERT,
     1                 TRIANG,TETRA,
     2                 RMAX,BASE,DEPTH,ZRATIO,X,Y,Z,XC,YC,ZC)
C
      IMPLICIT  NONE
      INTEGER   I,J,KK,NTET,INOD1,INOD2
      INTEGER   N,NT,NNOD,NSTR,NTRI,IPRT1,IVERT
      INTEGER   TRIANG(4,*),TETRA(5,*)
      REAL*8    ZMIN,ZTHICK,ZRSUM
      REAL*8    RMAX,BASE,DEPTH(*)
      REAL*8    ZRATIO(*),X(*),Y(*),Z(*),XC(*),YC(*),ZC(*)
      INCLUDE  'IOUNITS.H'
C
      open(99,file='output/grid3d')
      open(98,file='output/grid2d.exp')
      DO J=1,NSTR
         DO I=1,NTRI
            NTET=3*((J-1)*NTRI+I-1)
            TETRA(1,NTET+1)=(J-1)*NNOD+TRIANG(1,I)
            TETRA(2,NTET+1)=(J-1)*NNOD+TRIANG(2,I)
            TETRA(3,NTET+1)=(J-1)*NNOD+TRIANG(3,I)
            TETRA(4,NTET+1)=J*NNOD+TRIANG(1,I)
            TETRA(5,NTET+1)=TRIANG(4,I)
            TETRA(1,NTET+2)=TETRA(4,NTET+1)
            TETRA(2,NTET+2)=J*NNOD+TRIANG(2,I)
            TETRA(3,NTET+2)=J*NNOD+TRIANG(3,I)
            TETRA(4,NTET+2)=TETRA(3,NTET+1)
            TETRA(5,NTET+2)=TRIANG(4,I)
            TETRA(1,NTET+3)=TETRA(2,NTET+1)
            TETRA(2,NTET+3)=TETRA(3,NTET+1)
            TETRA(3,NTET+3)=TETRA(2,NTET+2)
            TETRA(4,NTET+3)=TETRA(4,NTET+1)
            TETRA(5,NTET+3)=TRIANG(4,I)
         END DO
      END DO
      ZMIN=RMAX
      DO I=1,NNOD
         IF (Z(I) .LT. ZMIN) ZMIN=Z(I)
      END DO
      DO I=1,NNOD
         ZTHICK=(Z(I) - ZMIN) + BASE
         ZRSUM=0.0D0
         DO J=1,NSTR
            KK=J*NNOD+I
            X(KK)=X(I)
            Y(KK)=Y(I)
            ZRSUM=ZRSUM + ZRATIO(J)
            IF (IVERT .EQ. 0) THEN
               Z(KK)=Z(I) - ZRSUM*BASE
            ELSE
               IF (IVERT .EQ. 1) THEN
                  Z(KK)=Z(I) - ZRSUM*ZTHICK
               ELSE IF (IVERT .EQ. 2) THEN
                  Z(KK)=ZMIN - ZRSUM*BASE
               ELSE IF (IVERT .EQ. 3) THEN
                  Z(KK)=Z(I) - ZRSUM*(DEPTH(I))
	       ELSE IF (IVERT .EQ. 4) THEN
	          Z(KK)=Z(I) - ZRSUM*BASE
	          IF (J .EQ. NSTR) THEN
	              Z(KK)=ZMIN-BASE
	          END IF
               END IF
            END IF
         END DO
      END DO
C     DO I=1,NDIR
C        CONTP(3,I)=CONTP2(I)
C     END DO
C     DO I=1,NSTR
C        DO J=1,NDIR
C           CONTP(3,I*NDIR+J)=NNOD*I+CONTP2(J)
C        END DO
C     END DO
      CALL NODELT(NT,TETRA,X,XC)
      CALL NODELT(NT,TETRA,Y,YC)
      CALL NODELT(NT,TETRA,Z,ZC)
      IF (IPRT1 .EQ. 3) THEN
         WRITE(IOUT2,1000)
         WRITE(ITERM,1000)
         WRITE(IOUT3,1010) NNOD,N
         DO I=1,N
            WRITE(IOUT3,1020) I,X(I),Y(I),Z(I)
         END DO
         WRITE(99,1030) NNOD,N,NT
         DO I=1,NT
            WRITE(99,1040) (TETRA(j,i),j=1,5)
         END DO
         DO I=1,N
            WRITE(99,1050) X(I),Y(I),Z(I)
         END DO
         WRITE(98,1060) NTRI,NNOD
         DO I=1,NNOD
            WRITE(98,1070) I,X(I),Y(I)
         END DO
         DO I=1,NTRI
            WRITE(98,1080) I,(TRIANG(j,i),j=1,4)
         END DO
         CALL CLOSIO
         STOP
      END IF
      close(99)
      close(98)
C
      RETURN
 1000 FORMAT(//,' IPRT1=3: Program terminating after output of X, Y, Z',
     1          ' coordinate values')
 1010 FORMAT(2I7,'  NNOD   N')
 1020 FORMAT(I7,3(1PE15.6))
 1030 FORMAT(3I9)
 1040 FORMAT(5I7)
 1050 FORMAT(3(1PE15.6))
 1060 FORMAT(2I7,1x,'1 0')
 1070 FORMAT('N ',I7,2(1PE15.6))
 1080 FORMAT('E ',5I7,'.')
 1100 FORMAT(//,' INPUT ERROR : ELEVATION VALUES NOT IN DESCENDING',
     1          ' ORDER ON SEEPAGE FACE ',I6,
     2       /,4X,'NODES',I4,' (NODE #',I6,') AND',I4,' (NODE #',I6,')')
      END
