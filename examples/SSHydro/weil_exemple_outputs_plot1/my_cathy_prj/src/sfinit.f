C
C**************************  SFINIT ************************************
C
C  calculate seepage face exit points based on initial pressure heads
C
C***********************************************************************
C
      SUBROUTINE SFINIT(N,NSF,NSFNUM,NSFNOD,SFEX,SFEXP,SFEXIT,DUPUIT,
     1                  PTIMEP,PNEW,SFFLAG,DELTAT,Z,PUNTDIRSF_NODE)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I,J,K,INOD,N
      INTEGER  NSF,DUPUIT
      INTEGER  NSFNUM(*),NSFNOD(NSFMAX,*)
      INTEGER  SFEX(NSFMAX,*),SFEXP(NSFMAX,*),SFEXIT(NSFMAX,*)
      INTEGER  SFFLAG(*)
      INTEGER  PUNTDIRSF_NODE(*) 
      REAL*8   DELTAT
      REAL*8   PTIMEP(*),PNEW(*),Z(*)
      INCLUDE 'IOUNITS.H'
C
      DO I=1,N
         PUNTDIRSF_NODE(I)=0 
      END DO
   
      WRITE(IOUT10,*) 'TIME=0.0d0'
      DO I=1,NSF
         DO J=1,NSFNUM(I)
            INOD=NSFNOD(I,J)
            IF (PTIMEP(INOD) .GE. 1.0e-8) THEN
               SFEX(I,J)=1
               IF (DUPUIT.EQ.0) THEN
                    PTIMEP(INOD)=0.0D0
                    PNEW(INOD)=0.0D0
               ELSE
                    WRITE(*,*)'ONLY DUPUIT = 0 FOR THIS CODE'
                    STOP
c                    PTIMEP(INOD)=0.0D0+Z(NSFNOD(I,SFEX(I)))-Z(INOD)
c                    PNEW(INOD)=0.0D0+Z(NSFNOD(I,SFEX(I)))-Z(INOD)
               END IF
C                   WRITE(666,*) INOD,PNEW(INOD)
            ELSE
               SFEX(I,J)=0
            END IF
            SFEXP(I,J)=SFEX(I,J)
            SFEXIT(I,J)=SFEX(I,J)
         END DO
         WRITE(IOUT10,*) (SFEX(I,J),J=1,NSFNUM(I))
      END DO
C     
      DO I=1,NSF
         DO J=1,NSFNUM(I)
            IF (SFEX(I,J) .EQ. 1) PUNTDIRSF_NODE(NSFNOD(I,J))=1
         END DO 
      END DO 
      RETURN
      END
