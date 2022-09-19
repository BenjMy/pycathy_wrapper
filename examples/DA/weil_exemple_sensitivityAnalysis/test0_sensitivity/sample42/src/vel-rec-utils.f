cccccc INIZIALIZZAZIONE VARIBILI GEOMETRICHE RELATIVE A
cccccc VELOCITY RECONSTRUCTION
c      SUBROUTINE CONNECTINITAL(TotalMeshFaces,Node_ElementCount,
c     1       Node_ElementIDs,Node_FaceCount,Node_FaceIDs,
c     1       Faces_ElementIDs,Faces_ElementCount,Faces_NodeIDs,
c     1       Faces_ReferenceVector,FaceCentroid,FaceArea,
c     1       Elements_FaceCount,Elements_FaceIDs,
c     1       ElementCentroid)

c     IMPLICIT NONE 

c      INCLUDE 'CATHY.H'
cc LOCALI
c      INTEGER I,J
cc Variabili da passare 
c      INTEGER  TotalMeshFaces,Node_ElementCount(NMAX)
c      INTEGER  Node_ElementIDs(NMAX,N1MAX),Faces_NodeIDs(NFACEMAX,3)
c      INTEGER  Faces_ElementCount(NFACEMAX),Faces_ElementIDs(NFACEMAX,2)
c      INTEGER  Node_FaceCount(NMAX),Node_FaceIDs(NMAX,MaxFConToNode)
c      INTEGER  Elements_FaceCount(NTEMAX),Elements_FaceIDs(NTEMAX,4)
c      REAL*8  Faces_ReferenceVector(NFACEMAX,3),FaceCentroid(NFACEMAX,3)
c      REAL*8   ElementCentroid(NTEMAX,3),FaceArea(NFACEMAX)
c      TotalMeshFaces=0
c
c      DO I=1,NMAX
c         Node_ElementCount(I)=0
c      END DO 

c     DO I=1,NMAX
c          DO J=1,N1MAX
c             Node_ElementIDs(I,J)=0
c          END DO 
c      END DO 
c
c      DO I=1,NFACEMAX
c         Faces_ElementCount(I)=0
c         DO J=1,2
c            Faces_ElementIds(I,J)=0
c         END DO 
c      END DO
c
c      DO I=1,NFACEMAX
c         DO J=1,3
c           Faces_NodeIDs(I,J)=0
c         END DO 
c      END DO 
c
c      DO I=1,NMAX
c         Node_FaceCount(I)=0
c      END DO
c
c      DO I=1,NMAX
c         DO J=1,MaxFConToNode
c            Node_FaceIDs(I,J)=0
c         END DO 
c      END DO   
       
c      DO I=1,NFACEMAX
c         DO J=1,3
c            Faces_ReferenceVector(I,J)=0.0d0
c         END DO 
c      END DO

c      DO I=1,NTEMAX
c         Elements_FaceCount(I)=0
c      END DO 

c      DO I=1,NTEMAX
c         DO J=1,4
c         Elements_FaceIDs(I,J)=0
c         END DO 
c      END DO

c      DO I=1,NFACEMAX
c         FaceArea(I)=0.0d0
c         DO J=1,3
c           FaceCentroid(I,J)=0.0d0
c         END DO 
c      END DO 
c
c      DO I=1,NTEMAX
c         DO J=1,3
c           ElementCentroid(I,J)=0.0d0
c         END DO 
c      END DO 
c
c      RETURN

c      END  


cc
      SUBROUTINE GetTet4FaceNodes(FaceNumber,Node1,Node2,Node3)

      IMPLICIT NONE
      
      INTEGER FaceNumber
      INTEGER Node1,Node2,Node3

      IF(FaceNumber.EQ.1) THEN
        Node1 = 1
        Node2 = 2
        Node3 = 3
      ELSE IF(FaceNumber.Eq.2) THEN
        Node1 = 1
        Node2 = 3
        Node3 = 4
      ELSE IF(FaceNumber.Eq.3) THEN
        Node1 = 1
        Node2 = 2
        Node3 = 4
      ELSE
        Node1 = 2
        Node2 = 3
        Node3 = 4
      END IF
      RETURN

      END SUBROUTINE

cccc  Crea la lista delle facce  
c      SUBROUTINE InsertFaceInList(FaceNode1,FaceNode2,FaceNode3,
c     1           Index,Element,Faces_ElementCount,Faces_NodeIDs,
c     1           Faces_ElementIDs,
c     1           Node_FaceCount,Node_FaceIDs,
c     1           TotalMeshFaces)
c
c      IMPLICIT NONE
c      INCLUDE 'CATHY.H'
c
c      INTEGER FaceNode1,FaceNode2,FaceNode3,Element
c      INTEGER Index
cc   Local varibles
c      LOGICAL Found
c      INTEGER Pos
cc  New variables
c      INTEGER TotalMeshFaces
c      INTEGER Faces_NodeIDs(NFACEMAX,3) 
c      INTEGER Node_FaceCount(NMAX)
c      INTEGER Node_FaceIDs(NMAX,MaxFConToNode)
c      INTEGER Faces_ElementCount(NFACEMAX)
c      INTEGER Faces_ElementIDs(NFACEMAX,2)
cc     INTEGER MaxFaces 
c
cc  mi serve per NMAX 
c
c      Found = .FALSE.
c      Index = 0
c
c      DO WHILE ((.NOT.Found).AND.(Index.LT.TotalMeshFaces))
cc       Update
c            Index = Index +1
c            Found = ((Faces_NodeIDs(Index,1).EQ.FaceNode1).AND.
c     1              (Faces_NodeIDs(Index,2).EQ.FaceNode2).AND.
c     1              (Faces_NodeIDs(Index,3).EQ.FaceNode3))
c      END DO
cc     Se non lo hai trovato aggiungilo
c      write(1111,*)FaceNode1,FaceNode2,FaceNode3
c      IF (.NOT.Found) THEN
c         IF ((TotalMeshFaces+1).LE.NFACEMAX) THEN
c            TotalMeshFaces =TotalMeshFaces + 1
c            Faces_NodeIDs(TotalMeshFaces,1) = FaceNode1
c            Faces_NodeIDs(TotalMeshFaces,2) = FaceNode2
c            Faces_NodeIDs(TotalMeshFaces,3) = FaceNode3
ccc    Restituisci l'indice della nuova faccia
c            Index = TotalMeshFaces
cc    Aggiungi la faccia alla lista delle facce afferenti al nodo
cc         Node 1
c            Node_FaceCount(FaceNode1) =
c     1      Node_FaceCount(FaceNode1) + 1
c            Pos = Node_FaceCount(FaceNode1)
c            Node_FaceIDs(FaceNode1,Pos) = TotalMeshFaces 
cc         Node 2
c            Node_FaceCount(FaceNode2) =
c     1      Node_FaceCount(FaceNode2) + 1
c            Pos = Node_FaceCount(FaceNode2)
c            Node_FaceIDs(FaceNode2,Pos) = TotalMeshFaces 
ccc         Node 3
c            Node_FaceCount(FaceNode3) =
c     1      Node_FaceCount(FaceNode3) + 1
c            Pos = Node_FaceCount(FaceNode3)
c            Node_FaceIDs(FaceNode3,Pos) = TotalMeshFaces 
ccc  aggiungi l'elemento alla lista di quelli afferenti alla faccia
c            Faces_ElementCount(TotalMeshFaces) = 
c     1      Faces_ElementCount(TotalMeshFaces) + 1 
c            Pos = Faces_ElementCount(TotalMeshFaces)
c            Faces_ElementIDs(TotalMEshFaces,Pos) = Element
c         ELSE
c           WRITE(*,*) 'Errore. Troppi spigoli. Programma Terminato'
c           STOP
c         END IF
c      ELSE
cc  se é gia' stato inseritoaggiungi solo l'elemento
c        Faces_ElementCount(Index) = 
c     1  Faces_ElementCount(Index) + 1 
c        Pos = Faces_ElementCount(Index)
c        Faces_ElementIDs(Index,Pos) = Element
c      END IF 
c      END SUBROUTINE

cccc CALCOLO del prodotto vettoriale tra due vettori
cccc (Per esempio per il vettore normale alla faccia)

      SUBROUTINE Do3DVecProduct(V1,V2,VP)

      IMPLICIT NONE 

      REAL*8 V1(3),V2(3),VP(3)

      VP(1) = V1(2)*V2(3)-V2(2)*V1(3)
      VP(2) = V1(3)*V2(1)-V2(3)*V1(1)
      VP(3) = V1(1)*V2(2)-V2(1)*V1(2)

      END SUBROUTINE 

cc CALCOLo prodotto scalare tra due vettori 
      REAL*8 FUNCTION Do3DScalarProduct(v,w)
      IMPLICIT NONE
      REAL*8 v(3),w(3)
      Do3DScalarProduct = v(1)*w(1)+v(2)*w(2)+v(3)*w(3)
      RETURN 
      END FUNCTION
cccccccccccc



cccc Per Normalizzare Un vettore
      SUBROUTINE Normalize3DVector(Versor)

      IMPLICIT NONE 

      REAL*8 Versor(3)
      REAL*8 Modulus,MathZero
cc MathZero da stabilire  
      MathZero = 1.0d-07
      Modulus = SQRT(Versor(1) * Versor(1) +
     r               Versor(2) * Versor(2) +
     r               Versor(3) * Versor(3))
      IF (ABS(Modulus).GT.MathZero) THEN
        Versor(1) = Versor(1)/Modulus
        Versor(2) = Versor(2)/Modulus
        Versor(3) = Versor(3)/Modulus
      ELSE
        Versor(1) = 0.0d0
        Versor(2) = 0.0d0
        Versor(3) = 0.0d0
      END IF
      END SUBROUTINE

CCC INIZIALIZZAZIONE VARIABILI RELTIVE ALLA 
CCC VELOCITY RECONSTRUCTION

      SUBROUTINE VELOCITYINITAL(NFACE,NT,UU,VV,WW,
     1           SIDE_CNC,PLIST,
     1           Faces_FaceType,
     1           ISIDE,
     1           XC,YC,ZC,FaceCentroid, 
     1           Faces_ReferenceVector, 
     1           Faces_NeumannFlux,
     1           Faces_FaceFlux)


      IMPLICIT NONE 
      INCLUDE 'CATHY.H'
cc LOCALI 
      INTEGER I,J,K,L
      INTEGER Node
      INTEGER CurrentFaceID
      INTEGER Faces_ElementCount(NFACEMAX) 
      REAL*8  CurrentNeumannFlux,CurrentFaceNormal(3)
      REAL*8  Flux(3)
cc Variabili Geometriche e relative ai flussi passate
cc ELEMENTI
      INTEGER NT
      INTEGER SIDE_CNC(4,NTEMAX)
      REAL*8  XC(NTEMAX),YC(NTEMAX),ZC(NTEMAX)
      REAL*8  FaceCentroid(3,NFACEMAX)
      REAL*8  UU(NTEMAX),VV(NTEMAX),WW(NTEMAX)
cc FACCE 
      INTEGER NFACE,Faces_FaceType(NFACEMAX)
      INTEGER ISIDE(3,NFACEMAX),PLIST(2,NFACEMAX)
      REAL*8  Faces_ReferenceVector(3,NFACEMAX)
ccc NODI
cc New variables da passare 
      REAL*8  Faces_FaceFlux(3,NFACEMAX),Faces_NeumannFlux(NFACEMAX)

cc Inizialize the 3 components of the flux for each face

      DO I=1,NFACE
         DO J=1,3
            Faces_FaceFlux(J,I) = 0.0d0
         END DO
      END DO
cc Assign Inital Elements Flux to each face
      DO I=1,NT
         Flux(1) = UU(I)
         Flux(2) = VV(I)
         Flux(3) = WW(I)
         DO J=1,4
cc Richiama la faccia globale
            CurrentFaceID=SIDE_CNC(J,I)
            IF((PLIST(2,CurrentFaceID).NE.0).
     1        OR.((PLIST(2,CurrentFaceID).EQ.0).AND.
     1        (Faces_FaceType(CurrentFaceID).EQ.0)))THEN
                DO K=1,3
                   Faces_FaceFlux(K,CurrentFaceID) = 
     1             Faces_FaceFlux(K,CurrentFaceID) + Flux(K)
                END DO
ccc CASO IN CUI LA FACCIA É DI BORDO   
            ELSE
            IF((PLIST(2,CurrentFaceID).EQ.0).AND.
     1         Faces_FaceType(CurrentFaceID).EQ.2)THEN
               CurrentNeumannFlux=Faces_NeumannFlux(CurrentFaceID)
            CALL EvalElementFaceNormal(I,J,CurrentFaceNormal,
     1           Faces_ReferenceVector,XC,YC,ZC,
     1           SIDE_CNC,FaceCentroid)
            DO K=1,3
                 Faces_FaceFlux(K,CurrentFaceID)=
     1           Faces_FaceFlux(K,CurrentFaceID)+
     1           CurrentNeumannFlux*CurrentFaceNormal(K)
            END DO
            END IF 
           END IF
         END DO 
       END DO  
c  Ciclo sui Nodi 
c               DO L=1,3
c                  Node=Faces_NodeIDs(CurrentFaceID,L)
ccc Nodo di Dirichelt 
c                  IF(Node_Type(Node).EQ.0)THEN
c                     DO K=1,3
c                        Faces_FaceFlux(CurrentFaceID,K) = 
c     1                  Faces_FaceFlux(CurrentFaceID,K)+Flux(K)/3.0d0
c                     END DO
c                  END IF
ccc Nodo di Neumann nullo 
c                  IF(Node_Type(Node).EQ.1)THEN
c                     DO K=1,3
c                       Faces_FaceFlux(CurrentFaceID,K) = 
c     1                 Faces_FaceFlux(CurrentFaceID,K) + 0.0d0
c                     END DO 
c                  END IF
cccc Neumann non Nullo. Qui c'é da pensare come portare i flussi
cc  volumetrici in flussi   
c                  IF(Node_Type(Node).EQ.2)THEN
c             Faces_Flux(CurrentFaceID)
c                     CurrentNeumannFlux = Node_NeumannFlux(Node)
c                     CALL EvalElementFaceNormal(I,J,CurrentFaceNormal,
c     1                    Faces_ReferenceVector,ElementCentroid,
c     1                    Elements_FaceIDs,FaceCentroid)
c                     DO K=1,3
c                        Faces_FaceFlux(CurrentFaceID,K)=
c     1                  CurrentNeumannFlux*
c     1                  CurrentFaceNormal(K)/3.0d0
c                     END DO
c                  END IF  
c               END DO 
c            END IF
ccc Divide by the number of elements sharing that edge 
       DO I=1,NFACE
          IF(PLIST(2,I).EQ.0)THEN
             Faces_ElementCount(I)=1
          ELSE
             Faces_ElementCount(I)=2
          END IF
          DO J=1,3
             Faces_FaceFlux(J,I) = 
     1       Faces_FaceFlux(J,I)/Faces_ElementCount(I)
          END DO 
       END DO 
      
       RETURN 
      
       END 


ccccc Per calcolare la normale alla faccia 
     
      SUBROUTINE EvalElementFaceNormal(CurrentElement,LocalFaceNumber,
     1           Normal,Faces_ReferenceVector,XC,YC,ZC,
     1           SIDE_CNC,FaceCentroid)

      IMPLICIT NONE
     
      INCLUDE 'CATHY.H'
cc Locali 
      INTEGER I,CurrentFaceID
      REAL*8  Normal(3),Aux(3)
      REAL*8  Do3DScalarProduct
      REAL*8  ElementCentroid(3,NTEMAX)
cc Variabili passate 
      INTEGER CurrentElement,LocalFaceNumber
      INTEGER SIDE_CNC(4,NTEMAX)
      REAL*8  XC(NTEMAX),YC(NTEMAX),ZC(NTEMAX),FaceCentroid(3,NFACEMAX) 
      REAL*8  Faces_ReferenceVector(3,NFACEMAX)
      REAL*8  Sign

      CurrentfaceID=SIDE_CNC(LocalFaceNumber,CurrentElement)
      Normal(1)=Faces_ReferenceVector(1,CurrentFaceID)
      Normal(2)=Faces_ReferenceVector(2,CurrentFaceID)
      Normal(3)=Faces_ReferenceVector(3,CurrentFaceID)
     
      ElementCentroid(1,CurrentElement)=XC(CurrentElement) 
      ElementCentroid(2,CurrentElement)=YC(CurrentElement) 
      ElementCentroid(3,CurrentElement)=ZC(CurrentElement) 
      DO I=1,3
         Aux(I)=FaceCentroid(I,CurrentFaceID)-
     1          ElementCentroid(I,CurrentElement)
      END DO

      CALL Normalize3DVector(Aux)
cccc DEvo considerare le normali esterne per questo faccio 
cccc questa operazione
      Sign = Do3DScalarProduct(Normal,Aux) 

      IF(Sign.LT.10.E-07)THEN
        DO I=1,3
           Normal(I)=-1.0d0*Normal(I)
c Non posso perché la CurrentFaceNormal é diversa quando 
c considero un elemento o l'altro 
c           Faces_ReferenceVector(CurrentFaceID,I)=Normal(I)
        END DO 
      END IF

      END SUBROUTINE 

CC  CARLOTTA 
      SUBROUTINE EvalElementFaceNormal2(CurrentElement,CurrentFaceID,
     1           Normal,Faces_ReferenceVector,XC,YC,ZC,
     1           FaceCentroid)

      IMPLICIT NONE
     
      INCLUDE 'CATHY.H'
cc Locali 
      INTEGER I,CurrentFaceID
      REAL*8  Normal(3),Aux(3)
      REAL*8  Do3DScalarProduct
      REAL*8  ElementCentroid(3,NTEMAX)
cc Variabili passate 
      INTEGER CurrentElement,LocalFaceNumber
      INTEGER SIDE_CNC(4,NTEMAX)
      REAL*8  XC(NTEMAX),YC(NTEMAX),ZC(NTEMAX),FaceCentroid(3,NFACEMAX) 
      REAL*8  Faces_ReferenceVector(3,NFACEMAX)
      REAL*8  Sign

c      CurrentfaceID=SIDE_CNC(LocalFaceNumber,CurrentElement)
      Normal(1)=Faces_ReferenceVector(1,CurrentFaceID)
      Normal(2)=Faces_ReferenceVector(2,CurrentFaceID)
      Normal(3)=Faces_ReferenceVector(3,CurrentFaceID)
     
      ElementCentroid(1,CurrentElement)=XC(CurrentElement) 
      ElementCentroid(2,CurrentElement)=YC(CurrentElement) 
      ElementCentroid(3,CurrentElement)=ZC(CurrentElement) 
      DO I=1,3
         Aux(I)=FaceCentroid(I,CurrentFaceID)-
     1          ElementCentroid(I,CurrentElement)
      END DO

      CALL Normalize3DVector(Aux)
cccc DEvo considerare le normali esterne per questo faccio 
cccc questa operazione
      Sign = Do3DScalarProduct(Normal,Aux) 

      IF(Sign.LT.10.E-07)THEN
        DO I=1,3
           Normal(I)=-1.0d0*Normal(I)
c Non posso perché la CurrentFaceNormal é diversa quando 
c considero un elemento o l'altro 
c           Faces_ReferenceVector(CurrentFaceID,I)=Normal(I)
        END DO 
      END IF

      END SUBROUTINE 

ccc  FUNCTION per il calcolo dell'imbalance su ciascun elemento
      REAL*8 FUNCTION EvalNetElementFlux(Element,PLIST, 
     1          SIDE_CNC,AREFACE,Faces_ReferenceVector,
     1          XC,YC,ZC,FaceCentroid,Faces_FaceFlux,
     1          UU,VV,WW)

      IMPLICIT NONE 

      INCLUDE 'CATHY.H'
ccc Variabili Locali 
      INTEGER I,J
      INTEGER CurrentFaceID
      REAL*8  NetElementFlux
      REAL*8  CurrentFaceFlux(3),CurrentFaceNormal(3)
      REAL*8  NormalFaceFlux
ccc FUNCTION 
      REAL*8  Do3DScalarProduct
ccc Variabili Passate 
      INTEGER Element
      INTEGER SIDE_CNC(4,NTEMAX)
      INTEGER PLIST(2,NFACEMAX)
      REAL*8  AREFACE(NFACEMAX),Faces_ReferenceVector(3,NFACEMAX)
      REAl*8  XC(NTEMAX),YC(NTEMAX),ZC(NTEMAX),FaceCentroid(3,NFACEMAX)
      REAL*8  Faces_FaceFlux(3,NFACEMAX)
      REAL*8  UU(NTEMAX),VV(NTEMAX),WW(NTEMAX) 

cccc Inizializzazione per ciascun elemento va inizializzata
      NetElementFlux=0.0d0 
      
      DO I=1,4
         CurrentFaceID=SIDE_CNC(I,Element)
ccc Considero solo le facce interne perché globalmente il flusso 
ccc dovrebbe annullarsi
         IF(PLIST(2,CurrentFaceID).NE.0)THEN
               CALL EvalElementFaceNormal(Element,I,CurrentFaceNormal,
     1              Faces_ReferenceVector,XC,YC,ZC,
     1              SIDE_CNC,FaceCentroid)
         
            CurrentFaceFlux(1)=UU(Element)
            CurrentFaceFlux(2)=VV(Element)
            CurrentFaceFlux(3)=WW(Element)
         
         NormalFaceFLux=Do3DScalarProduct(CurrentFaceFlux,
     1             CurrentFaceNormal)*AREFACE(CurrentFaceID) 

         NetElementFlux = NetElementFlux + NormalFaceFLux
         END IF
      END DO
      EvalNetElementFlux = NetElementFlux 
    
      END FUNCTION 

ccc Make a copy of the Faces_FaceFlux before post-processing
      
c--------------dopo la soluzione ricalcolo dei nuovi flussi----------- 
!     
! UPDATE FACE FLUXES WITH CORRECTIONS
!----------------------------------------------------------------------
      SUBROUTINE UpdateFaceFluxes(NFACE,Faces_FaceFlux,
     1           Faces_FluxCorrection,Faces_ReferenceVector,
     1           PLIST,
     1           N,Node_FaceCount,Node_FaceIDs,
     1           SIDE_CNC,
     1           XC,YC,ZC,FaceCentroid,UU,VV,WW)
      IMPLICIT NONE
!     Include Global Variables
      INCLUDE 'CATHY.H'
      INTEGER NFACE,N
      INTEGER Node_FaceCount(NMAX)
      INTEGER Node_FaceIDs(NMAX,N1MAX)
      INTEGER PLIST(2,NFACEMAX)
      INTEGER SIDE_CNC(4,NTEMAX)
      REAL*8  Faces_ReferenceVector(3,NFACEMAX)
      REAL*8  XC(NTEMAX),YC(NTEMAX),ZC(NTEMAX),FaceCentroid(3,NFACEMAX)
      REAl*8  Faces_FaceFlux(3,NFACEMAX)
      REAl*8  Faces_FluxCorrection(NFACEMAX,3)
      REAl*8  UU(NTEMAX),VV(NTEMAX),WW(NTEMAX)
!     Local Variables
      INTEGER I,J,K
      INTEGER CurrentFaceID,CurrentElement
      INTEGER LocalFaceNumber
      REAl*8  CurrentNeumannFlux,Normal(3),Flux(3)
!     Routine Core
      DO I = 1,NFACE
!       Update
cccc Faces_FLuxCorrection sono tre scalari, li moltiplico per il
ccc Faces_RefVector e ottengo il vettore
        DO J = 1,3
           Faces_FaceFlux(J,I) =
     1     Faces_FaceFlux(J,I) + (1.0d0/3.0d0)*
     1     (Faces_FluxCorrection(I,1)+
     1     Faces_FluxCorrection(I,2)+
     1     Faces_FluxCorrection(I,3))*
     1     Faces_ReferenceVector(J,I)
        END DO
!       Reset Corrections
        DO J = 1,3
          Faces_FluxCorrection(I,J) = 0.0d0
        END DO
      END DO
c      DO I=1,N
c        IF(Node_Type(I).NE.3)THEN 
c          DO J=1,Node_FaceCount(I)
c           CurrentFaceID=Node_FaceIDs(I,J)
c           IF(Faces_ElementCount(CurrentFaceID).EQ.1)THEN
c              CurrentElement=Faces_ElementIDs(CurrentFaceID,1)
cccc Trovo il LocalFaceNumber
c             DO K=1,4
c               IF(Elements_FaceIDs(CurrentElement,K).
c     1            EQ.CurrentFaceID)THEN
c                 LocalFaceNumber=K  
c               END IF   
c             END DO  
c       CALL EvalElementFaceNormal(CurrentElement,LocalFaceNumber,Normal,
c     1          Faces_ReferenceVector,ElementCentroid,
c     1          Elements_FaceIDs,FaceCentroid)
cccc  Ora devi sommare i flussi 
c            IF(Node_Type(I).EQ.0)THEN 
c             Flux(1)=UU(CurrentElement)
c              Flux(2)=VV(CurrentElement)
c              Flux(3)=WW(CurrentElement)
c              DO K=1,3
c                 Faces_FaceFlux(CurrentFaceID,K)=
c     1           Faces_FaceFlux(CurrentFaceID,K)+Flux(K)
c              END DO
c            ELSE IF(Node_Type(I).EQ.2)THEN
c              CurrentNeumannFlux=Node_NeumannFlux(I)
c              DO K=1,3
c                 Faces_FaceFlux(CurrentFaceID,K)=
c     1           Faces_FaceFLux(CurrentFaceID,K)+CurrentNeumannFlux*
c     1           Normal(K)
c              END DO 
c            END IF  
c           END IF
c         END DO 
c        END IF
c       END DO  
      RETURN
      END SUBROUTINE

cc-----------------------------------------------------------------
      SUBROUTINE EvalCentroidElementFluxes(NT, 
     1           AREFACE,XC,YC,ZC,VOLU,TETRA,
     1           SIDE_CNC,
     1           Elements_FluxXCentroid,Elements_FluxYCentroid,
     1           Elements_FluxZCentroid,Faces_ReferenceVector,
     1           FaceCentroid,X,Y,Z,ISIDE,
     1           Faces_FaceFlux)
      IMPLICIT NONE
!     Include Global Variables
      INCLUDE 'CATHY.H'
c   variabili da passare 
      INTEGER NT
      INTEGER TETRA(5,NTEMAX),ISIDE(3,NFACEMAX)
      INTEGER SIDE_CNC(4,NTEMAX)
      REAL*8  AREFACE(NFACEMAX),XC(NTEMAX),YC(NTEMAX),ZC(NTEMAX)
      REAl*8  VOLU(NTEMAX),Elements_FluxXCentroid(NTEMAX)
      REAL*8  Elements_FLuxYCentroid(NTEMAX)
      REAl*8  Elements_FLuxZCentroid(NTEMAX)
      REAL*8  Faces_ReferenceVector(3,NFACEMAX)
      REAl*8  FaceCentroid(3,NFACEMAX)
      REAL*8  X(NMAX),Y(NMAX),Z(NMAX)
      REAl*8  Faces_FaceFlux(3,NFACEMAX)
!     Local Variables
      INTEGER I
!     Sub Core
      DO I = 1,NT
        CALL Update3DCentroidElementFlux(I,AREFACE,XC,YC,ZC,
     1       VOLU,Elements_FLuxXCentroid,Elements_FluxYCentroid,
     1       Elements_FluxZCentroid,
     1       SIDE_CNC,
     1       Faces_ReferenceVector,FaceCentroid,X,Y,Z,TETRA,
     1       ISIDE,Faces_FaceFlux)
      END DO
c      write(201,*) Elements_FluxXCentroid(50),
c     1             Elements_FluxYCentroid(50),
c     1             Elements_FluxZCentroid(50)
      
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE Update3DCentroidElementFLux(CurrentElement,AREFACE,
     1           XC,YC,ZC,VOLU,Elements_FluxXCentroid,
     1           Elements_FluxYCentroid,Elements_FluxZCentroid,
     1           SIDE_CNC,
     1           Faces_ReferenceVector,
     1           FaceCentroid,X,Y,Z,TETRA,ISIDE,
     1           Faces_FaceFlux)
      IMPLICIT NONE
!     Include Global Variables
      INCLUDE 'CATHY.H'
!     Arguments
      INTEGER CurrentElement
!     External Functions
      INTEGER SIDE_CNC(4,NTEMAX)
      INTEGER TETRA(5,NTEMAX),ISIDE(3,NFACEMAX)
      REAL*8  Faces_ReferenceVector(3,NFACEMAX)
      REAL*8  AREFACE(NFACEMAX),XC(NTEMAX),YC(NTEMAX),ZC(NTEMAX)
      REAL*8  FaceCentroid(3,NFACEMAX)
      REAL*8  VOLU(NTEMAX),X(NMAX),Y(NMAX),Z(NMAX)
      REAL*8  Elements_FLuxXCentroid(NTEMAX)
      REAl*8  Elements_FluxYCentroid(NTEMAX)
      REAl*8  Elements_FluxZCentroid(NTEMAX)
      INTEGER EvalNodeOppositeToFace
      REAL*8  Do3DScalarProduct
      REAL*8  Faces_FaceFlux(3,NFACEMAX)
!     Local Variables
      INTEGER J,K
      REAL*8  CurrentElementCentroid(3)
      REAL*8  CurrentElementVolume,CurrentFaceArea
      INTEGER CurrentFaceID
      REAL*8  CurrentFaceNormal(3)
      INTEGER CurrentNode
      REAL*8  CurrentVector(3)
      REAL*8  CurrentFaceFlux(3)
      REAL*8  CurrentFaceNormalFlux
      REAL*8  RT0_ShapeValue(3)
!     Routine Core
!      Initialize
       Elements_FluxXCentroid(CurrentElement) = 0.0d0
       Elements_FluxYCentroid(CurrentElement) = 0.0d0
       Elements_FluxZCentroid(CurrentElement) = 0.0d0
!      Eval Area
       CurrentElementVolume = VOLU(CurrentElement)
          CurrentElementCentroid(1)=XC(CurrentElement)
          CurrentElementCentroid(2)=YC(CurrentElement)
          CurrentElementCentroid(3)=ZC(CurrentElement)
       DO J = 1,4
        CurrentFaceID = SIDE_CNC(J,CurrentElement)
!       Valuta la lunghezza dello spigolo
        CurrentFaceArea = AREFACE(CurrentFaceID)
!       Eval Normal
        CALL EvalElementFaceNormal(CurrentElement,J,CurrentFaceNormal,
     1           Faces_ReferenceVector,XC,YC,ZC,
     1           SIDE_CNC,FaceCentroid)
c        CALL EvalElementFaceNormal(CurrentElement,J,
c     r                             CurrentFaceNormal)
!       Trova il Nodo Opposto allo spigolo
        CurrentNode = EvalNodeOppositeToFace(CurrentFaceID,
     1                CurrentElement,TETRA,ISIDE)
!       Eval Distance
        CurrentVector(1) = CurrentElementCentroid(1)-
     1                     X(CurrentNode)
        CurrentVector(2) = CurrentElementCentroid(2)-
     1                     Y(CurrentNode)
        CurrentVector(3) = CurrentElementCentroid(3)-
     1                     Z(CurrentNode)
c  FINE
!       Eval Current Edge Flux
        DO K = 1,3
          CurrentFaceFlux(K) =
     1    Faces_FaceFlux(K,CurrentFaceID)
        END DO
!       Eval Normal Flux
        CurrentFaceNormalFlux = 
     1  Do3DScalarProduct(CurrentFaceFlux,
     1                    CurrentFaceNormal)*CurrentFaceArea
!       Eval Raviart-Thomas Shape Value
        DO K = 1,3
           RT0_ShapeValue(K) = (1.0d0/(3.0d0*CurrentElementVolume))*
     1                            CurrentVector(K)
        END DO
!       Sum up Values
c        DO LoopC = 1,3
          Elements_FluxXCentroid(CurrentElement) =
     1    Elements_FluxXCentroid(CurrentElement) +
     1    RT0_ShapeValue(1)*CurrentFaceNormalFlux
          Elements_FluxYCentroid(CurrentElement) =
     1    Elements_FluxYCentroid(CurrentElement) +
     1    RT0_ShapeValue(2)*CurrentFaceNormalFlux
          Elements_FluxZCentroid(CurrentElement) =
     1    Elements_FluxZCentroid(CurrentElement) +
     1    RT0_ShapeValue(3)*CurrentFaceNormalFlux
c        END DO
       END DO
c  Quebsta non l'ho capita
       DO K = 1,3
          CurrentFaceFlux(K) = Faces_FaceFlux(K,CurrentFaceID)
       END DO
      RETURN 

      END SUBROUTINE
c-----------------------------------------------
      INTEGER  FUNCTION EvalNodeOppositeToFace(CurrentFaceID,
     1                   CurrentElement,TETRA,ISIDE)
      IMPLICIT NONE
!     Include Global Variables
      INCLUDE 'CATHY.H'
ccc variabili passate e da passare 
      INTEGER ISIDE(3,NFACEMAX) 
!     Function Arguments
      INTEGER CurrentFaceID,CurrentElement
      INTEGER TETRA(5,NTEMAX)
!     Local Variables
      INTEGER*4 Count
      LOGICAL Found
      INTEGER*4 LocalResult
!     Function Core
      Count = 0
      Found = .FALSE.
      DO WHILE ((.NOT.Found).AND.(Count.LT.4))
        Count = Count + 1
        LocalResult = TETRA(Count,CurrentElement)
        Found =((LocalResult.NE.ISIDE(1,CurrentFaceID)).AND.
     r          (LocalResult.NE.ISIDE(2,CurrentFaceID)).AND.
     r          (LocalResult.NE.ISIDE(3,CurrentFaceID)))
      END DO
      IF (.NOT.Found) THEN
        WRITE(*,*)'Problema Interno in EvalNodeOppositeToFace.'
        STOP
      END IF
      EvalNodeOppositeToFace = LocalResult
      RETURN
      END FUNCTION
!-----------------------------------------------------------------------
