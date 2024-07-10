C
C**************************  EXTCVG ************************************
C
C  check for convergence of seepage face exit points 
C
C***********************************************************************
C
      SUBROUTINE EXTCVG(NSF,NSFNOD,SFEX,SFEXIT,TIME,DELTAT,ITER,KSF,
     1                  NSFNUM)
C
      IMPLICIT  NONE
      INCLUDE  'CATHY.H'
      INTEGER   I,J
      INTEGER   NSF,ITER,KSF
      INTEGER   NSFNUM(*)
      INTEGER   NSFNOD(NSFMAX,*),SFEX(NSFMAX,*),SFEXIT(NSFMAX,*)
      REAL*8    TIME,DELTAT
      INCLUDE  'IOUNITS.H'
C
      KSF=0
      DO I=1,NSF
         DO J=1,NSFNUM(I)
            IF (SFEX(I,J) .NE. SFEXIT(I,J)) KSF=KSF + 1
         END DO
      END DO
c      WRITE(IOUT10,1100) TIME,DELTAT,ITER
c      WRITE(IOUT10,1105) 
c      WRITE(IOUT10,1115) (SFEXIT(I),NSFNOD(I,SFEXIT(I)),I=1,NSF)
c      WRITE(IOUT10,1110) 
c      WRITE(IOUT10,1115) (SFEX(I),NSFNOD(I,SFEX(I)),I=1,NSF)
      IF (KSF .EQ. 0) THEN
c         WRITE(IOUT10,1120) 
      ELSE
c         WRITE(IOUT10,1125)
      END IF
C
      RETURN
 1100 FORMAT(/,' (TIME=',1PE12.4,', DELTAT=',1PE12.4,', ITER=',I5,')')
 1105 FORMAT(  1X,'EXIT POINTS PREV ITER: NODE (NODE #)')
 1110 FORMAT(  1X,'EXIT POINTS CURR ITER: NODE (NODE #)')
 1115 FORMAT(6(2X,I3,'(',I6,')'))
 1120 FORMAT(1X,'SEEPAGE FACE EXIT POINTS CONVERGED',/)
 1125 FORMAT(1X,'SEEPAGE FACE EXIT POINTS DID NOT CONVERGE',/)
      END
