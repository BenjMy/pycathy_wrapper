CCC   With this subroutine I am calculating the type of face 
CCC   necessary for velocity reconstruction. 
CCC   IF Faces_FaceType=-1 internal face
CCC   IF Faces_FaceType=2  Neumann Face
CCC   IF Faces_FaceType=0  Dirichlet Face (fixed head)
CCC   Note*--It is a Dirichlet face if it has one of the three nodes 
CCC   of Dirichlet type. 
CCC   (In reality the type of face should 
CCC   be set depending on which local node we are considering 
CCC   in the reconstruction: if we are performing the reconstruction 
CCC   around a Dirichlet node, the boundary faces connected to that 
CCC   node are of Dirichlet type, of Neumann otherwise--to be extended) 
CCC   Note**-Atmospheric nodes are always considered of Neumann
CCC   Sign Convention for Faces_NeumannFlux-Negative if of inflow
CCC
      SUBROUTINE UPD_FACETYPE(N,NNOD,NFACE,PLIST,ISIDE,Faces_FaceType,
     1                       NODE_FACECOUNT,NODE_FACEIDS,
     2                       ATMACT,ARENOD,AREFACE,Faces_NeumannFLux,
     3                       PUNTDIRFLOW_NODE,IFSF,PUNTDIRSF_NODE,IFATM,
     4                       PUNTNEUFLOW_NODE,CONNEU_NODE)

       IMPLICIT NONE
       INCLUDE 'CATHY.H'
cc LOCAL
       INTEGER I,K
       INTEGER Node1,Node2,Node3,Faccia
       REAL*8  AREA(NMAX)
cc INPUT 
       INTEGER N,NNOD,NFACE
       INTEGER PUNTDIRFLOW_NODE(*)
       INTEGER PUNTDIRSF_NODE(*)
       INTEGER NODE_FACECOUNT(*),NODE_FACEIDS(NMAX,*) 
       INTEGER IFSF(*)
       INTEGER PLIST(2,*),ISIDE(3,*)
       INTEGER IFATM(*),PUNTNEUFLOW_NODE(*) 
c
       REAL*8 ARENOD(*)
       REAL*8 ATMACT(*),AREFACE(*)
       REAL*8 CONNEU_NODE(*) 
c  OUTPUT 
       INTEGER Faces_FaceType(*)
       REAL*8  Faces_NeumannFLux(*) 
C
C   First I am setting all faces internal
       DO I=1,NFACE
          Node1=ISIDE(1,I)
          Node2=ISIDE(2,I)
          Node3=ISIDE(3,I)
          AREA(Node1)=0.0d0
          AREA(Node2)=0.0d0
          AREA(Node3)=0.0d0
          Faces_NeumannFlux(I)=0.0d0
C  Internal Face 
          Faces_FaceType(I)=-1
C  If They Are Boundary Faces 
          IF(PLIST(2,I).EQ.0)THEN
            Faces_FaceType(I)=2
C  Neumann Face (Neumann imposed conditions are of flow L^3/T). 
C  For this reason it is necessary to weight them on the area. 
C  (it is the flux L/T that can be interpolated)
            DO K=1,Node_FaceCount(Node1)
               Faccia=Node_FaceIDS(Node1,K)
               IF(PLIST(2,Faccia).EQ.0)THEN 
                 AREA(Node1)=AREA(Node1)+AREFACE(Faccia)/3.0d0
               END IF
            END DO 
            DO K=1,Node_FaceCount(Node2)
               Faccia=Node_FaceIDS(Node2,K)
               IF(PLIST(2,Faccia).EQ.0)THEN 
                 AREA(Node2)=AREA(Node2)+AREFACE(Faccia)/3.0d0
               END IF
            END DO 
            DO K=1,Node_FaceCount(Node3)
               Faccia=Node_FaceIDS(Node3,K)
               IF(PLIST(2,Faccia).EQ.0)THEN 
                 AREA(Node3)=AREA(Node3)+AREFACE(Faccia)/3.0d0
               END IF
            END DO 
            Faces_NeumannFlux(I)=-(CONNEU_NODE(Node1)/AREA(Node1)+
     1      CONNEU_NODE(Node2)/AREA(Node2)+
     2      CONNEU_NODE(Node3)/AREA(Node3))/3.0d0
C  Dirichlet Face if it has a Dirichlet node 
            IF((PUNTDIRFLOW_NODE(Node1).EQ.1).OR.
     1        (PUNTDIRFLOW_NODE(Node2).EQ.1).OR.
     2        (PUNTDIRFLOW_NODE(Node3).EQ.1))THEN
               Faces_FaceType(I)=0
            END IF
C Seepage face faces if they have all three nodes of SF              
            IF((IFSF(Node1).EQ.1).AND.
     1         (IFSF(Node2).EQ.1).AND.
     2         (IFSF(Node3).EQ.1))THEN
C Dirichlet Face if one of the three nodes is of Dirichlet, 
C otherwise it stays of 0-flux Neumann
                IF((PUNTDIRSF_NODE(Node1).EQ.1).OR.
     1             (PUNTDIRSF_NODE(Node2).EQ.1).OR.
     2             (PUNTDIRSF_NODE(Node3).EQ.1))THEN
                   Faces_FaceType(I)=0
                END IF
            END IF
          END IF 
       END DO
C       
C Imposing Neumann condition on atmospheric nodes. If IFATM of one of 
C the surface node is atmospheric (IFATM.NE.-1) then all surface 
C nodes are atmospheric 
C       
       IF(IFATM(1).NE.-1)THEN 
         DO I=1,NFACE
            Node1=ISIDE(1,I)
            Node2=ISIDE(2,I)
            Node3=ISIDE(3,I)
            IF((Node1.LE.NNOD).AND.
     1        (Node2.LE.NNOD).AND.
     2        (Node3.LE.NNOD))THEN
              Faces_NeumannFlux(I)=-(ATMACT(Node1)/ARENOD(Node1)+
     1        ATMACT(Node2)/ARENOD(Node2)+
     1        ATMACT(Node3)/ARENOD(Node3))/3.0d0
*AREFACE(I)
            END IF
         END DO  
       END IF

       RETURN 

       END 
