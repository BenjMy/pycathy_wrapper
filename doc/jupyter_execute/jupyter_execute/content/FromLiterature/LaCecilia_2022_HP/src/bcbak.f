C
C**************************  BCBAK *************************************
C
C  special handling of time variable non-atmospheric, non-seepage face
C  boundary conditions for the case where BCTIM(1) < TIME <= BCTIM(2)
C  during back-stepping 
C
C***********************************************************************
C
      SUBROUTINE BCBAK(BCTYPE,IOUNI2,
     1                 IPRT1,NBC,TIME,BCTIM,BCINP,NNBC,BC,ANBC,ANNBC)
C
      IMPLICIT     NONE
      INTEGER      I,ANBC
      INTEGER      IOUNI2,IPRT1
      INTEGER      NBC(3)
      INTEGER      ANNBC(*)
      INTEGER      NNBC(3,*)
      REAL*8       SLOPE,TIM21R,TIMW1
      REAL*8       TIME
      REAL*8       BCTIM(*),BCINP(3,*),BC(*)
      CHARACTER*22 BCTYPE
C
C  we don't need to do anything if first time step or if
C  non-atmospheric, non-seepage face boundary conditions inputs are 
C  homogeneous in time
C
      IF (BCTIM(1) .GE. BCTIM(2)) GO TO 800
C     TIM21R=1.0D0/(BCTIM(2) - BCTIM(1))
C     TIMW1=TIME - BCTIM(1)
      ANBC=NBC(1)
      DO I=1,ANBC
C        SLOPE=(BCINP(2,I) - BCINP(1,I))*TIM21R
C        BC(I)=BCINP(1,I) + SLOPE*TIMW1
         BC(I)=BCINP(1,I)
	 ANNBC(I)=NNBC(1,I)
      END DO
      IF (IPRT1 .GE. 1) THEN
         WRITE(IOUNI2,1060) BCTYPE,TIME
         WRITE(IOUNI2,1010) (I,ANNBC(I),BC(I),I=1,ANBC)
      END IF
C
  800 RETURN
 1010 FORMAT((1X,3(I5,I6,1X,1PE13.5)))
 1060 FORMAT(2X,A22,' BCs VALUES AT TIME: ',1PE12.4,/,
     1       1X,3('    #  NODE  BC VALUE    '))
      END
