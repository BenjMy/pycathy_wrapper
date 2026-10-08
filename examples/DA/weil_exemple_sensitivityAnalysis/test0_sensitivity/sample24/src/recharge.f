C
C**************************  RECHARGE **********************************
C
C  calculates water table recharge in terms of both flux and volume
C
C***********************************************************************
C
      SUBROUTINE RECHARGE(NNOD,NSTR,WNOD,ARENOD,PSI,TIME,DELTAT,RECFLOW,
     1                    RECVOL,RECNOD)
C
      IMPLICIT  NONE
      INTEGER   I,J
      INTEGER   NNOD,NSTR
      REAL*8    TIME,DELTAT,ZERO
      REAL*8    RECFLOW,RECVOL
      REAL*8    WNOD(*),ARENOD(*),PSI(*)
      REAL*8    RECNOD(*)
      PARAMETER (ZERO=0.0d0)
C
      RECFLOW = zero
      IF (TIME.EQ.0.0d0) THEN
         RECVOL = zero
      END IF
      CALL INIT0R(NNOD,RECNOD)
C
      DO I=NNOD*NSTR+1,NNOD*(NSTR+1)
         DO J=1,NSTR
            IF (PSI(I-(J-1)*NNOD).GT.ZERO.AND.PSI(I-J*NNOD).LE.ZERO) 
     1      THEN
               IF (WNOD(I-J*NNOD).LE.ZERO) THEN
                  RECNOD(I-NNOD*NSTR)=-1.0d0*WNOD(I-J*NNOD)*
     1                                  ARENOD(I-NSTR*NNOD)
                  GO TO 100
               END IF
            END IF
         END DO
         IF (PSI(I-NNOD*NSTR).GE.0.0d0) THEN
            IF (WNOD(I-NNOD*NSTR).LE.ZERO) THEN
               RECNOD(I-NNOD*NSTR)=-1.0d0*WNOD(I-NSTR*NNOD)*
     1                              ARENOD(I-NSTR*NNOD)
            END IF
         END IF
 100     RECFLOW=RECFLOW+RECNOD(I-NNOD*NSTR)
      END DO
      RECVOL=RECVOL+RECFLOW*DELTAT
      return
      end
