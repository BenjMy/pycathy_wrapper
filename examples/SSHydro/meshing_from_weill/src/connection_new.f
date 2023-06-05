C-----------------------------------------------------------------------
      SUBROUTINE CONNECTION_NEW(N,NT,NFACE,TETRA,SIDE_CNC,PLIST,ISIDE,
     1     NEIGH,AREFACE,FaceCentroid,Node_ElementCount,Node_ElementIDs,
     2     NODE_FACECOUNT,NODE_FACEIDS,
     3     Faces_ReferenceVector,X,Y,Z)
      
C-----------------------------------------------------------------------
C  Subroutine for computing the connection between nodes and faces
C
C  Input Variables:
C
C  ntet
C  ntetmax
C  nfacmax
C  tetra(4,ntetmax)
C
C  Output Variables:
C
C  nfac
C  side_cnc(4,ntet)
C  neigh(4,ntet)
C  plist(2,nfac)
C  iside(3,nfac)
C  
C  Work arrays:
C
C  faces(5,4*ntetmax)
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      include 'CATHY.H'
C
C  Input Variables
C
      integer  nt,n
      integer  tetra(5,ntemax)
      real*8   x(*),y(*),z(*)
C     
C  Output Variables
C
      integer  nface
      integer  side_cnc(4,ntemax)
      integer  neigh(4,ntemax)
      integer  plist(2,nfacemax),iside(3,nfacemax)
      integer  node_elementcount(nmax),node_elementids(nmax,n1max)
      integer  node_facecount(nmax),node_faceids(nmax,MaxFConToNode)
      real*8   areface(nfacemax)
      real*8   Faces_ReferenceVector(3,*)
      real*8   FaceCentroid(3,*)
C
C  Work arrays
C
      integer  faces(5,4*ntemax)
      integer  point(nfacemax)

C
C  Local variables
C
      integer  perm(3,4)
      integer  i,j,k,jj,kk
      integer  imax
      integer  pos,DumInt
      integer  Node1,Node2,Node3
      integer  FaceNode1,FaceNode2,FaceNode3
c  Local Real variables
      real*8   x1,x2,x3,y1,y2,y3,z1,z2,z3
      real*8   s1,s2,s3
      real*8   Ref(3),Aux1(3),Aux2(3)



      character*100 file_topol

      data imax/2147483647/
      data perm/1,2,3,1,2,4,1,3,4,2,3,4/

C
C  Generate faces vector

      do i = 1,nt
C        Sort the element nodes in increasing order
         call sort_elem(4,tetra(1,i))
         do j = 1,4
            kk = 4*(i-1) + j
            faces(1,kk) = tetra(perm(1,j),i)
            faces(2,kk) = tetra(perm(2,j),i)
            faces(3,kk) = tetra(perm(3,j),i)
            faces(4,kk) = i
            faces(5,kk) = 0
         enddo
      enddo
C
C  Sort in lexicographic order the faces vector
C
      call heapsort(faces,4*nt)

c
c  Eliminate doubled faces and build side_cnc
c
C     Set side_cnc and neigh to 0
      do j=1,nt
         do i=1,4 
            side_cnc(i,j)=0
            neigh(i,j)=0
         end do
      end do
c
      nface=1
      do i=2,4*nt
c
c        Check the faces, if coincident remove the newer
c        and update faces(5,nface) so that faces(4,nface) < faces(5,nface)
c
         if ((faces(1,i).eq.faces(1,nface)).and. 
     1       (faces(2,i).eq.faces(2,nface)).and. 
     2       (faces(3,i).eq.faces(3,nface))) then
            if (faces(4,i).gt.faces(4,nface)) then
               faces(5,nface) = faces(4,i)
            else 
               faces(5,nface) = faces(4,nface)
               faces(4,nface) = faces(4,i)
            end if
         else
           nface=nface+1
           do j=1,5
              faces(j,nface)=faces(j,i)
           end do
         endif
      end do
c
c  Build the neigh vector using plist as a work array
c
      do i = 1,nface
         point(i) = 1
      enddo

      do i = 1,nface
         iside(1,i)=faces(1,i)
         iside(2,i)=faces(2,i)
         iside(3,i)=faces(3,i)
         plist(1,i)=faces(4,i)
         plist(2,i)=faces(5,i)
         jj = faces(4,i)
         kk = faces(5,i)
         if (kk .ne. 0) then
            neigh(point(jj),jj) = kk
            neigh(point(kk),kk) = jj
            point(jj) = point(jj) + 1
            point(kk) = point(kk) + 1            
         endif
      enddo
c
      do i=1,nface
         j=faces(4,i)
         k=faces(5,i)
         kk=1
         do while (kk.le.4)
            if (side_cnc(kk,j).eq.0) then
               side_cnc(kk,j)= i
               kk=5
            end if
            kk=kk+1
         end do
         if (k.gt.0) then
            kk=1
            do while (kk.le.4)              
               if (side_cnc(kk,k).eq.0) then
                  side_cnc(kk,k)= i
                  kk=5
               end if
               kk=kk+1
            end do
         end if
      end do
      DO I=1,NFACE
           x1=x(iside(1,I))
           x2=x(iside(2,I))
           x3=x(iside(3,I))
           y1=y(iside(1,I))
           y2=y(iside(2,I))
           y3=y(iside(3,I))
           z1=z(iside(1,I))
           z2=z(iside(2,I))
           z3=z(iside(3,I))
           s1=(x1*(y2-y3)-y1*(x2-x3)+(x2*y3-y2*x3))
           s2=(x1*(z2-z3)-z1*(x2-x3)+(x2*z3-z2*x3))
           s3=(y1*(z2-z3)-z1*(y2-y3)+(y2*z3-z2*y3))
           AREFACE(I)=0.5d0*DSQRT(s1*s1+s2*s2+s3*s3)
      END DO 

C COMPUTE THE NODE_ELEMENT CONNECTIVITIES
        DO I=1,N
           Node_ElementCount(I)=0
           DO J=1,N1MAX
              Node_ElementIDs(I,J)=0
           END DO  
        END DO 
        DO I=1,NT
           DO J=1,4
              CALL GetTet4FaceNodes(J,Node1,Node2,Node3)
              FaceNode1 = TETRA(Node1,I)
              FaceNode2 = TETRA(Node2,I)
              FaceNode3 = TETRA(Node3,I)
c Aggiungi l'elemento corrente alla lista di quelli afferenti al nodo 
              Node_ElementCount(TETRA(J,I))=
     1        Node_ElementCount(TETRA(J,I)) + 1
cc Assegna la numerazione globale all'elemento numerato localmente
              Pos = Node_ElementCount(TETRA(J,I)) 
              Node_ElementIDs(TETRA(J,I),Pos) = I
cc Ordina i nodi delle facce 
              IF (FaceNode1.GT.FaceNode2) THEN
                 DumInt = FaceNode2
                 FaceNode2 = FaceNode1
                 FaceNode1 = DumInt
              END IF
              IF (FaceNode1.GT.FaceNode3) THEN
                 DumInt = FaceNode3
                 FaceNode3 = FaceNode1  
                 FaceNode1 = DumInt
              END IF
              IF (FaceNode2.GT.FaceNode3) THEN
                 DumInt = FaceNode2
                 FaceNode2 = FaceNode3
                 FaceNode3 = DumInt
              END IF          
c             CALL InsertFaceInList(FaceNode1,FaceNode2,FaceNode3,
c     1            FaceListIndex,I,Faces_ElementCount,Faces_NodeIDs,
c     1            Faces_ElementIDs,
c     1            Node_FaceCount,Node_FaceIDs,
c     1            TotalMeshFaces) 
!  Aggiungi l'indice alla lista delle facce per elementi
c                  Elements_FaceCount(I) =
c     1            Elements_FaceCount(I) + 1
c                  Pos =Elements_FaceCount(I)
c                  Elements_FaceIDs(I,Pos)= FaceListIndex
           END DO 
        END DO 
C CARLOTTA Inizialization 

      DO I=1,N
         Node_FaceCount(I)=0
         DO J=1,MaxFConToNode
            Node_FaceIDs(I,J)=0
         END DO 
      END DO       
      DO I=1,NFACE
         DO J=1,3
            Node_FaceCount(ISIDE(J,I)) =
     1      Node_FaceCount(ISIDE(J,I)) + 1
            Pos = Node_FaceCount(ISIDE(J,I))
            Node_FaceIDs(ISIDE(J,I),Pos) = I 
         END DO 
      END DO 
C   Faces_ReferenceVector (required for VELOCITY RECONSTRUCTION)
        DO I=1,NFACE
           DO J=1,3
              Faces_ReferenceVector(J,I)=0.0d0
           END DO 
        END DO 
        DO I = 1,NFACE 
           Node1 = ISIDE(1,I)     
           Node2 = ISIDE(2,I)     
           Node3 = ISIDE(3,I) 
c  Aux1
           Aux1(1) = X(Node2) - X(Node1)
           Aux1(2) = Y(Node2) - Y(Node1)
           Aux1(3) = Z(Node2) - Z(Node1)
c  Aux2  
           Aux2(1) = X(Node3) - X(Node1)
           Aux2(2) = Y(Node3) - Y(Node1)
           Aux2(3) = Z(Node3) - Z(Node1)
c Calcolo la normale alla faccia come prodotto 
c vettoriale tra i due vettori
           CALL Do3DVecProduct(Aux1,Aux2,Ref) 
         
           CALL Normalize3DVector(Ref) 

           Faces_ReferenceVector(1,I) = Ref(1)
           Faces_ReferenceVector(2,I) = Ref(2)
           Faces_ReferenceVector(3,I) = Ref(3)
         END DO
C  Compute Face Centroid
         DO I=1,NFACE
            DO J=1,3
               FaceCentroid(J,I)=0.0d0
            END DO 
         END DO 
         DO I=1,NFACE
            DO J=1,3
               FaceCentroid(1,I)=FaceCentroid(1,I)+
     1         X(ISIDE(J,I)) 
               FaceCentroid(2,I)=FaceCentroid(2,I)+
     1         Y(ISIDE(J,I)) 
               FaceCentroid(3,I)=FaceCentroid(3,I)+
     1         Z(ISIDE(J,I)) 
            END DO
               FaceCentroid(1,I)=FaceCentroid(1,I)/3.0d0 
               FaceCentroid(2,I)=FaceCentroid(2,I)/3.0d0 
               FaceCentroid(3,I)=FaceCentroid(3,I)/3.0d0 
         END DO
       
      return
      end
c      SUBROUTINE GetTet4FaceNodes(FaceNumber,Node1,Node2,Node3)
c
c      IMPLICIT NONE
c      
c      INTEGER FaceNumber
c      INTEGER Node1,Node2,Node3
c
c      IF(FaceNumber.EQ.1) THEN
c        Node1 = 1
c        Node2 = 2
c        Node3 = 3
c      ELSE IF(FaceNumber.Eq.2) THEN
c        Node1 = 1
c        Node2 = 3
c        Node3 = 4
c      ELSE IF(FaceNumber.Eq.3) THEN
c        Node1 = 1
c        Node2 = 2
c        Node3 = 4
c      ELSE
c        Node1 = 2
c        Node2 = 3
c        Node3 = 4
c      END IF
c      RETURN
c
c      END SUBROUTINE
