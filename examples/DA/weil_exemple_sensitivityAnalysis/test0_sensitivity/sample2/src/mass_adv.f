C
C**************************  MASS_ADV ************************************
C
C  post-calculate the mass advected in and out after the 
C  advective router (VOLFIN) 
C  
C***********************************************************************
C
      SUBROUTINE MASS_ADV(NFACE,DELTATADV,PLIST,VN,
     1           PRECL,PRECR,MASSOUTADV,MASSINADV,
     2           MASSOUTTOTADV,MASSINTOTADV)
C
      IMPLICIT NONE
      INTEGER  I
      INTEGER  NFACE,PLIST(2,*)
      REAL*8   VN(*),DELTATADV
      REAL*8   PRECL(*),PRECR(*) 
      REAL*8   MASSOUTADV,MASSINADV
      REAL*8   MASSOUTTOTADV,MASSINTOTADV
      REAL*8   MASSINLOC,MASSOUTLOC
C
C  CARLOTTA POST COMPUTATION OF MASSBALANCE 
c      SFMASSOUT=0.0d0
c      SFEV=0.0d0 
c      SFIN=0.0d0
      MASSINLOC=0.0d0
      MASSOUTLOC=0.0d0
      DO I=1,NFACE
         IF(PLIST(2,I).EQ.0)THEN
           IF(VN(I).GE.0.0d0)THEN
CC QUI É SU TUTTO BOUNDARY
              MASSOUTADV=MASSOUTADV+VN(I)*PRECL(I)*DELTATADV
              MASSOUTLOC=MASSOUTLOC+VN(I)*PRECL(I)*DELTATADV
c           ELSE 
c             SFEV=SFEV+VN(I)*PRECL(I)*
c     1                 DELTATADV
c           END IF
           ELSE
              MASSINADV=MASSINADV+VN(I)*PRECR(I)*DELTATADV
              MASSINLOC=MASSINLOC+VN(I)*PRECR(I)*DELTATADV
           END IF
         END IF  
      END DO
C Condition bidouille, mais sinon il y a un
C problème d'écriture dans les fichiers de sortie
C     IF (MASSOUTLOC.LE.1.0D-99) THEN 
C     MASSOUTLOC=0.0d0
C     END IF
C     IF (MASSINLOC.LE.1.0D-99) THEN 
C     MASSINLOC=0.0d0
C     END IF     
C
c      SFMASSOUTTOT=SFMASSOUTTOT+SFMASSOUT
      MASSOUTTOTADV=MASSOUTTOTADV+MASSOUTLOC
c      SFEVTOT=SFEVTOT+SFEV
      MASSINTOTADV=MASSINTOTADV+MASSINLOC
cc  CARLOTTA Queste sono di Dirichlet       
c         DO I=1,NFACE
c            IF(Faces_FaceType(I).EQ.0)THEN
c              IF(VN(I).GT.0.0d0)THEN
cc qui sto guardando Seepage face
c                C=C+VN(I)*precl(i)*deltatadv
c              END IF
c           END IF
c         END DO   
C
      RETURN
      END
