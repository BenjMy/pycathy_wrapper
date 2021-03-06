       SUBROUTINE DIRFACE_UPDATE(NNOD,NFACE,ISIDE,ATMINT,ARENOD,
     1            CONCNOD_PERDIR_INT,NPFA_TRA,CONTPFA_TRA,
     2            PRESCFA_TRA,PUNTDIR)

       IMPLICIT NONE 
       INCLUDE 'CATHY.H'
c Variables passed
       INTEGER NNOD,NFACE,ISIDE(3,*)
       REAL*8  ATMINT(*),ARENOD(*),CONCNOD_PERDIR_INT(*)
c Local variables
       INTEGER I
       INTEGER N1,N2,N3
       REAL*8  ATMACT_FACE(NFACEMAX)
c Output variables
       INTEGER NPFA_TRA,CONTPFA_TRA(*)
       INTEGER PUNTDIR(*)
       REAL*8  PRESCFA_TRA(*)

       CALL INIT0R(NFACEMAX,ATMACT_FACE)

       NPFA_TRA=0 
       DO I=1,NFACE
          N1=ISIDE(1,I) 
          N2=ISIDE(2,I) 
          N3=ISIDE(3,I)
          IF((N1.LE.NNOD).AND.(N2.LE.NNOD).
     1     AND.(N3.LE.NNOD))THEN
         ATMACT_FACE(I)=(ATMINT(N1)/ARENOD(N1)+ATMINT(N2)/ARENOD(N2)+
     1                    ATMINT(N3)/ARENOD(N3))/3.0d0
         IF(ATMACT_FACE(I).GT.0.0d0)THEN 
           NPFA_TRA=NPFA_TRA+1
           CONTPFA_TRA(NPFA_TRA)=I
           PRESCFA_TRA(NPFA_TRA)=(ATMINT(N1)*CONCNOD_PERDIR_INT(N1)/
     1     ARENOD(N1)+CONCNOD_PERDIR_INT(N2)*ATMINT(N2)/ARENOD(N2)
     2     +ATMINT(N3)*CONCNOD_PERDIR_INT(N3)/ARENOD(N3))/
     3     3.0d0/ATMACT_FACE(I)
          END IF
        END IF
      END DO

      CALL BOUNDIR(NFACE,NPFA_TRA,CONTPFA_TRA,PUNTDIR)   

      RETURN

      END      
