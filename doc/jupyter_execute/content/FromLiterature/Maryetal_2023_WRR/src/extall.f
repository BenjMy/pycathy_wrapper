C
C**************************  EXTALL ************************************
C
C  calculate new position of the exit point along each seepage face,
C  and adjust boundary conditions for the seepage face nodes to reflect
C  changes in the position of the exit point.
C  Check for convergence of seepage face exit points. 
C  Case ISFONE=0 : seepage face exit point updating performed by 
C                  checking all nodes on a seepage face
C
C***********************************************************************
C
      SUBROUTINE EXTALL(N,NSF,NSFNUM,NSFNOD,SFEX,SFEXIT,SFQ,PNEW,
     1                  SFFLAG,ITER,TIME,DELTAT,KSF,DUPUIT,Z,
     2                  PUNTDIRSF_NODE)
C
      IMPLICIT  NONE
      INCLUDE  'CATHY.H'
      INTEGER   I,J,K,KEX,INOD
      INTEGER   N,NSF,ITER,KSF,DUPUIT
      INTEGER   NSFNUM(*),NSFNOD(NSFMAX,*),SFEX(NSFMAX,*)
      INTEGER   SFEXIT(NSFMAX,*),SFFLAG(*),PUNTDIRSF_NODE(*)
      REAL*8    TIME,DELTAT
      REAL*8    SFQ(NSFMAX,*),PNEW(*),Z(*)
      INCLUDE  'IOUNITS.H'
C
      DO I=1,N
         PUNTDIRSF_NODE(I)=0
      END DO
      WRITE(IOUT10,*) 'TIME=',TIME
      DO I=1,NSF
         DO J=1,NSFNUM(I)
            INOD=NSFNOD(I,J)
            IF (SFEX(I,J).EQ.1) THEN
               IF (SFQ(I,J) .GE. 0.0D0) THEN
                  SFEX(I,J)=0
                  SFQ(I,J)=0.0  
               END IF
            ELSE
               IF (PNEW(INOD) .GT. 1.0E-8) THEN
                  PNEW(INOD)=0.0D0
                  SFEX(I,J)=1
               END IF
            END IF
         END DO
         WRITE(IOUT10,*) (SFEX(I,J),J=1,NSFNUM(I))
      END DO

c     IF (DUPUIT.EQ.1) THEN
c        DO I=1,NSF
c         DO J=SFEX(I),NSFNUM(I)
c            PNEW(INOD)=0.0D0+Z(NSFNOD(I,SFEX(I)))-Z(INOD)
c         END DO
c        END DO
c     END IF
C
      DO I=1,NSF
         DO J=1,NSFNUM(I)
            IF (SFEX(I,J).EQ.1) PUNTDIRSF_NODE(NSFNOD(I,J))=1
         END DO
      END DO
C
C  check for convergence of seepage face exit points 
C
      CALL EXTCVG(NSF,NSFNOD,SFEX,SFEXIT,TIME,DELTAT,ITER,KSF,NSFNUM)
C
C     WRITE(6666,*)
      RETURN
 2100 FORMAT(  ' SFFLAG(1) AT SEEPAGE FACE ',I4,
     1         ' (TIME=',1PE8.2,', DELTAT=',1PE8.2,', ITER=',I4,')',
     2       /,11X,'NO EXIT POINT; SEEPAGE FACE COMPLETELY',
     3         ' UNSATURATED')
 2200 FORMAT(  ' SFFLAG(2) AT SEEPAGE FACE ',I4,
     1         ' (TIME=',1PE8.2,', DELTAT=',1PE8.2,', ITER=',I4,')',
     2       /,11X,'EXIT POINT AT TOP OF SEEPAGE FACE; SEEPAGE FACE',
     3         ' COMPLETELY SATURATED')
 2300 FORMAT(  ' SFFLAG(3) AT SEEPAGE FACE ',I4,
     1         ' (TIME=',1PE8.2,', DELTAT=',1PE8.2,', ITER=',I4,')',
     2       /,11X,'NON-NEGATIVE PRESSURE HEAD OF ',1PE9.3,' AT NODE ',
     3         I3,' (NODE # ',I6,')',
     4       /,11X,'OCCURRED AT A POTENTIAL SEEPAGE FACE NODE DURING',
     5         ' EXIT POINT LOWERING')
 2400 FORMAT(  ' SFFLAG(4) AT SEEPAGE FACE ',I4,
     1         ' (TIME=',1PE8.2,', DELTAT=',1PE8.2,', ITER=',I4,')',
     2       /,11X,'NON-NEGATIVE PRESSURE HEAD OF ',1PE9.3,' AT NODE ',
     3         I3,' (NODE # ',I6,')',
     4       /,11X,'OCCURRED AT A POTENTIAL SEEPAGE FACE NODE DURING',
     5         ' EXIT POINT RAISING')
      END
