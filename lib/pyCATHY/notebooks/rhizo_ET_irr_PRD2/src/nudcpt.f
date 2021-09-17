C
C**************************  NUDCPT ************************************
C
C  evaluate model (computed) results at the observation points.
C  i.e., for those elements that contain observation points, interpolate
C  the computed soil moisture values on the nodes of these elements
C  to the coordinates of the observation point.
C
C***********************************************************************
C
      SUBROUTINE NUDCPT(NUDFLAG,NUDN,TETRA,NUDTET,IVOL,VOLUR,SW,PNODI,
     1                  PTNEW,AI,BI,CI,DI,NUDX,NUDY,NUDZ,NUDSMC)
C
      IMPLICIT NONE
      INTEGER  I,IN,IEL,INOD
      INTEGER  NUDN,NUDFLAG
      INTEGER  TETRA(5,*),NUDTET(*),IVOL(*)
      REAL*8   VOLUR(*),SW(*),PNODI(*),PTNEW(*)
      REAL*8   AI(4,*),BI(4,*),CI(4,*),DI(4,*)
      REAL*8   NUDX(*),NUDY(*),NUDZ(*),NUDSMC(*)
C
      DO I=1,NUDN
         IEL=NUDTET(I)
	 NUDSMC(I)=0.0D0
	 DO IN=1,4
	    INOD=TETRA(IN,IEL)
            IF (NUDFLAG .EQ. 0) THEN
               NUDSMC(I)=NUDSMC(I) +
     1                   (AI(IN,IEL) + BI(IN,IEL)*NUDX(I) +
     2                    CI(IN,IEL)*NUDY(I) + DI(IN,IEL)*NUDZ(I)) *
     3                   SW(INOD) * PNODI(INOD) 
            ELSE
               NUDSMC(I)=NUDSMC(I) +
     1                   (AI(IN,IEL) + BI(IN,IEL)*NUDX(I) +
     2                    CI(IN,IEL)*NUDY(I) + DI(IN,IEL)*NUDZ(I)) *
     3                   PTNEW(INOD)
            END IF
         END DO
         NUDSMC(I) = NUDSMC(I) * VOLUR(IEL) * IVOL(IEL)
      END DO
C
      RETURN
      END
