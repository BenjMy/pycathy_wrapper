C
C**************************  BCNXT_TRA  ************************************
C
C  input (if necessary) and interpolate time variable TRANSPORT DIRICHLET bc
C
C***********************************************************************
C
      SUBROUTINE BCNXT_TRA(BCTYPE,IUNIT,NNOD,NBCMAX,
     1                 IPRT1,NB2C,NBFC,NBC,NSTR,HTIBC,TIME,
     1                 BCTIM,BCINP,BCNOD,BC,ANBC,ABCNOD)
C
      IMPLICIT     NONE
      INTEGER      I,ANBC,NNOD,NBCMAX
      INTEGER      IUNIT,IPRT1,NSTR,HTIBC
      INTEGER      NB2C(3),NBFC(3),NBC(3) 
      INTEGER      ABCNOD(NBCMAX),BCNOD(3,NBCMAX)
      REAL*8       TIMEIN,SLOPE,TIM32R,TIMW1
      REAL*8       TIME
      REAL*8       BCTIM(3),BCINP(3,NBCMAX),BC(NBCMAX)
      CHARACTER*22 BCTYPE
C     INCLUDE 'CATHY.H'
C
      IF (HTIBC .NE. 0) GO TO 800
  200 IF (TIME .LE. BCTIM(3)) GO TO 300
      BCTIM(1)=BCTIM(2)
      BCTIM(2)=BCTIM(3)
      NB2C(1) = NB2C(2)
      NB2C(2) = NB2C(3)
      NBFC(1) = NBFC(2)
      NBFC(2) = NBFC(3)
      NBC(1) = NBC(2)
      NBC(2) = NBC(3)
      DO I=1,NBCMAX
         BCINP(1,I)=BCINP(2,I)
         BCINP(2,I)=BCINP(3,I)
	 BCNOD(1,I)=BCNOD(2,I)
	 BCNOD(2,I)=BCNOD(3,I)
      END DO
C      IF (IPRT1 .GE. 1) THEN
C         WRITE(IOUNI1,1050) BCTYPE,BCTIM(2)
C         WRITE(IOUNI2,1055) BCTYPE,BCTIM(2)
C         WRITE(IOUNI2,1010) (I,BCNOD(2,I),BCINP(2,I),I=1,NBC(2))
C      END IF
      READ(IUNIT,*,END=700) TIMEIN
      BCTIM(3)=TIMEIN
      CALL RDNDBC(IUNIT,NB2C,NBFC,NBC,NSTR,BCNOD,NNOD)
      CALL READBC(IUNIT,NB2C(3),NBFC(3),NBC(3),NSTR,NNOD,BCINP)
      GO TO 200
  300 IF (BCTIM(3) .GT. BCTIM(2)) THEN
C        TIM32R=1.0D0/(BCTIM(3) - BCTIM(2))
C        TIMW1=TIME - BCTIM(2)
         ANBC=NBC(2)
         DO I=1,ANBC
C           SLOPE=(BCINP(3,I) - BCINP(2,I))*TIM32R
C           BC(I)=BCINP(2,I) + SLOPE*TIMW1
            BC(I)=BCINP(2,I)
	    ABCNOD(I)=BCNOD(2,I)
         END DO
      ELSE
         ANBC=NBC(3)
         DO I=1,ANBC
            BC(I)=BCINP(3,I)
	    ABCNOD(I)=BCNOD(3,I)
         END DO
      END IF
C     WRITE(999,*) ANBC
C     WRITE(999,*)(ABCNOD(I),I=1,ANBC)
C     WRITE(999,*)(BC(I),I=1,ANBC)
C      IF (IPRT1 .GE. 1) THEN
C         WRITE(IOUNI2,1060) BCTYPE,TIME
C         WRITE(IOUNI2,1010) (I,BCNOD(3,I),BC(I),I=1,NBC(3))
C      END IF
      GO TO 800
C
  700 HTIBC=1
      GO TO 300
C
  800 RETURN
 1010 FORMAT((1X,3(I5,I6,1X,1PE13.5)))
 1050 FORMAT(/,2X,A22,' BCs INPUT AT TIME: ',1PE12.4)
 1055 FORMAT(2X,A22,' BCs  INPUT AT TIME: ',1PE12.4,/,
     1       1X,3('    #  NODE  BC VALUE    '))
 1060 FORMAT(2X,A22,' BCs VALUES AT TIME: ',1PE12.4,/,
     1       1X,3('    #  NODE  BC VALUE    '))
      END
