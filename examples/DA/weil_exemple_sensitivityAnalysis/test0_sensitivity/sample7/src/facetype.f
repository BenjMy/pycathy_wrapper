      SUBROUTINE FACETYPE(N,NFACE,PLIST,ISIDE,IFSF,IFSFP,IFATM,
     1                    IFATMP,IBNDRYP,IBNDRY,
     2                    NODE_FACECOUNT,NODE_FACEIDS,AREFACE)
C
      IMPLICIT NONE 
      INCLUDE 'CATHY.H' 
C  Local 
      INTEGER I
c  Input 
      INTEGER N,NFACE,PLIST(2,*),ISIDE(3,*)
      INTEGER IFSFP(*),IFATMP(*)
      INTEGER IFSF(*),IFATM(*)
      INTEGER NODE_FACEIDS(NMAX,*),NODE_FACECOUNT(*)
      REAL*8  AREFACE(*)
c  Output
      INTEGER IBNDRYP(*),IBNDRY(*)

      DO I=1,NFACE
         IF(PLIST(2,I).EQ.0)THEN
           IBNDRY(I)=0
         ELSE 
           IBNDRY(I)=-1
         END IF
      END DO
C SEARCH FOR THE SEEPAGE FACE FACES
      DO I=1,NFACE
         IF(PLIST(2,I).EQ.0)THEN
              IF((IFSF(ISIDE(1,I)).EQ.1).and.
     1          (IFSF(ISIDE(2,I)).EQ.1).AND.
     2          (IFSF(ISIDE(3,I)).EQ.1))THEN
                IBNDRY(I)=1      
              END IF
         END IF 
      END DO 
C AT THE PREVIOUS TIME LEVEL  
      DO I=1,NFACE
         IF(PLIST(2,I).EQ.0)THEN
              IF((IFSFP(ISIDE(1,I)).EQ.1).and.
     1          (IFSFP(ISIDE(2,I)).EQ.1).AND.
     2          (IFSFP(ISIDE(3,I)).EQ.1))THEN
                IBNDRYP(I)=1      
              END IF
         END IF 
      END DO
C  
c      DO I=1,N
c         ARENODO(I)=0.0d0
c      END DO  
c      DO I=1,N
c         IF(IFSFP(I).EQ.1)THEN
c           DO J=1,NODE_FACECOUNT(I)
c              IF(IBNDRYP(NODE_FACEIDS(I,J)).EQ.1)THEN
cc                 ARENODO(I)=ARENODO(I)+
c     1           AREFACE(NODE_FACEIDS(I,J))/3.0d0
c              END IF 
c           END DO 
c         END IF 
c      END DO 

      RETURN 

      END  

