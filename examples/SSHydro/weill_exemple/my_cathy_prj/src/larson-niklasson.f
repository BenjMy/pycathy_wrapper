ccc       larson-niklasson main 

          SUBROUTINE ApplyLN(N,NT,Node_ElementCount,Node_ElementIDs,
     1               TETRA,UU,VV,WW,X,Y,Z,VOLU,ET1E,PEL,PNEW,
     1               PTIMEP,SW_NODE,SWP_NODE,DELTAT,
     1               ISIDE,AREFACE,
     1               Faces_ReferenceVector,Faces_FaceFlux,
     1               XC,YC,ZC,FaceCentroid,SIDE_CNC,
     1               Faces_FaceType,Node_FaceIDs,Node_FaceCount,
     1               PLIST,
     1               Faces_FluxCorrection,
     1               Faces_NeumannFlux)

          IMPLICIT NONE 

          INCLUDE 'CATHY.H'
ccc locali 
          INTEGER I,J,K,L
          INTEGER CurrentFaceID
          INTEGER NodeIDx,CurrentElement
          REAl*8  CurrentSign,CurrentProduct
          REAL*8  ResidualVector(N1MAX)
          REAL*8  LocMat(N1MAX,N1MAX)
          REAl*8  LocalRef(3),Normal(3)
          REAl*8  MathZero
          REAl*8  Do3DScalarProduct
ccc passate 
          INTEGER N,NT,TETRA(5,NTEMAX)
          INTEGER Node_ElementCount(NMAX)
ccc passate e da passare 
          INTEGER Node_ElementIDs(NMAX,N1MAX)
          REAL*8  UU(NTEMAX),VV(NTEMAX),WW(NTEMAX)
          REAL*8  X(NMAX),Y(NMAX),Z(NMAX)
          REAL*8  XC(NTEMAX),YC(NTEMAX),ZC(NTEMAX)
          REAL*8  VOLU(NTEMAX)
ccc passate e da passare per il calcolo del residual
          REAL*8  Faces_NeumannFlux(NFACEMAX)
ccc passate e da passare per il calcolo di locmas
          INTEGER Faces_FaceType(NFACEMAX)
          INTEGER PLIST(2,NFACEMAX)
          INTEGER Node_FaceCount(NMAX)
          INTEGER Node_FaceIDs(NMAX,MaxFConToNode)
cc passate e da passare per il solutore 
cc passate e da passare relative al residual (derivata temporale)
          REAL*8  ET1E(NTEMAX),PEL(NTEMAX)
          REAL*8  PNEW(NMAX),PTIMEP(NMAX),SW_NODE(NMAX),SWP_NODE(NMAX)
          REAl*8  DELTAT
ccc Passate per il calcolo dei flusso non conservativo 
ccc assemblato con il residual
          INTEGER ISIDE(3,NFACEMAX)
          INTEGER SIDE_CNC(4,NTEMAX)
          REAL*8  AREFACE(NFACEMAX),Faces_ReferenceVector(3,NFACEMAX)
          REAL*8  Faces_FaceFlux(3,NFACEMAX) 
          REAL*8  FaceCentroid(3,NFACEMAX),ElementCentroid(3,NTEMAX)
cccc Soluzione del sistema
          REAl*8  Faces_FluxCorrection(NFACEMAX,3)
ccc Inizializzazione di LocMat e ReidualVector
         MathZero=1.0E-07
         DO I=1,NT
            Elementcentroid(1,I)=XC(I)
            Elementcentroid(2,I)=YC(I)
            Elementcentroid(3,I)=ZC(I)
         END DO
         DO I=1,N
             DO J=1,Node_ElementCount(I)
c             DO J=1,N1MAX
                ResidualVector(J)=0.0d0
                DO K=1,Node_ElementCount(I)
c                DO K=1,N1MAX
                   LocMat(J,K)=0.0d0
                END DO 
             END DO
ccc Assembla il termine noto per ciascun nodo

             CALL AssembleResidual(I,ResidualVector,Node_ElementCount,
     1            Node_ElementIDs,TETRA,UU,VV,WW,X,Y,Z,XC,YC,ZC,VOLU,
     1            ET1E,PEL,PNEW,PTIMEP,SW_NODE,SWP_NODE,DELTAT,
     1            ISIDE,
     1            AREFACE,Faces_ReferenceVector,
     1            Faces_FaceFlux,ElementCentroid,
     1            FaceCentroid,SIDE_CNC,
     1            PLIST)
c    Assembla la matrice dei coefficienti 
             CALL AssembleMatrix(I,LocMat,AREFACE,XC,YC,ZC,
     1            Node_ElementIDs,Node_ElementCount,
     1            SIDE_CNC,Node_FaceIDs,Node_FaceCount,
     1            PLIST,Faces_FaceType)
c    Risolvi il sistema locale 
!    Risolvi il Problema algebrico applicando le opportune condensazioni
             CALL SolveFor3DAuxiliaryFluxes(I,LocMat,ResidualVector,
     1            Node_ElementCount,Node_FaceCount,PLIST,
     1            Node_FaceIDs)
!    Salva le correzioni normali sulle facce
          DO J = 1,Node_ElementCount(I)
          CurrentElement = Node_ElementIDs(I,J)
           DO K = 1,4
            CurrentFaceID =
     1      SIDE_CNC(K,CurrentElement)
cccccc  Se la faccia non é di boundary 
            IF((PLIST(2,CurrentFaceID).NE.0).OR.
     1        ((PLIST(2,CurrentFaceID).EQ.0).AND.
     1        (Faces_FaceType(CurrentFaceID).EQ.0)))THEN
c            IF ((Faces_FaceType(CurrentFaceID).EQ.0).OR.
c     1         (Faces_FaceType(CurrentFaceID).EQ.3)) THEN
              IF ((ISIDE(1,CurrentFaceID).EQ.I).OR.
     1            (ISIDE(2,CurrentFaceID).EQ.I).OR.
     1            (ISIDE(3,CurrentFaceID).EQ.I)) THEN
c    Da qui 
!               Eval Edge Normal
            CALL EvalElementFaceNormal(CurrentElement,K,
     1           Normal,Faces_ReferenceVector,XC,YC,ZC,
     1           SIDE_CNC,FaceCentroid)
c            END IF
c                CALL EvalElementFaceNormal(CurrentElement,LoopC,
c     r                                     CurrentFaceNormal)
!               Fai un Copia Locale del Vettore di riferimento
                DO L = 1,3
                  LocalRef(L) =
     1            Faces_ReferenceVector(L,CurrentFaceID)
                END DO
!               Valuta il prodotto interno
cccc pensa bene a questa cosa perché non mi é chiara
                CurrentProduct =
     1          Do3DScalarProduct(Normal,LocalRef)
                IF (CurrentProduct.GT.MathZero)THEN
                  CurrentSign = 1.0d0
                ELSE
                  CurrentSign = -1.0d0
                END IF
!               Applica la Correzione
                IF (ISIDE(1,CurrentFaceID).EQ.I)THEN
                  NodeIdx = 1
                ELSE IF (ISIDE(2,CurrentFaceID).EQ.
     r                   I) THEN
                  NodeIdx = 2
                ELSE
                  NodeIdx = 3
                END IF
                Faces_FluxCorrection(CurrentFaceID,NodeIdx) =
     r          Faces_FluxCorrection(CurrentFaceID,NodeIdx) +
     r          CurrentSign*ResidualVector(J)
              END IF
              END IF
           END DO 
          END DO
         END DO
c          DO J=1,Node_ElementCount(I)
c             DO K=1,Node_ElementCount(I)
c                write(*,*)J,K,Locmat(J,K)
c              END DO 
c         END DO  
          RETURN 

          END

ccc--------------------ASSEMBLARESIDUAL--------------------------------
ccc 
        SUBROUTINE AssembleResidual(Node,ResidualVector,
     1             Node_ElementCount,Node_ElementIDs,TETRA,
     1             UU,VV,WW,X,Y,Z,XC,YC,ZC,VOLU,
     1             ET1E,PEL,PNEW,PTIMEP,SW_NODE,SWP_NODE,
     1             DELTAT,
     1             ISIDE,
     1             AREFACE,Faces_ReferenceVector,
     1             Faces_FaceFlux,ElementCentroid,
     1             FaceCentroid,SIDE_CNC,
     1             PLIST)
    
        IMPLICIT NONE 
        INCLUDE 'CATHY.H'
ccc Locali
        INTEGER J,K,L
        INTEGER CurrentElement,CurrentFaceID
        REAL*8  CurrentFaceArea,FaceNormalFlux
        REAL*8  CurrentNeumannFLux
        REAL*8  Normal(3),LocalFaceFlux(3)
c commento sono i nodi relativi all'elemento considerato, localnode number
c é invece il numero locale (1-4) del nodo dell'elemento considerato che 
c é = al nodo considerato
        INTEGER LocalNodeNumber 
ccc passate e da passare 
        REAL*8  UU(NTEMAX),VV(NTEMAX),WW(NTEMAX)
        REAL*8  VOLU(NTEMAX)
        REAL*8  X(NMAX),Y(NMAX),Z(NMAX)
        REAL*8  XC(NTEMAX),YC(NTEMAX),ZC(NTEMAX)
ccc Passate 
        INTEGER Node,Node_ElementCount(NMAX)
        INTEGER Node_ElementIDs(NMAX,N1MAX)
        INTEGER PLIST(2,NFACEMAX)
        INTEGER TETRA(5,NTEMAX)
        REAL*8  ResidualVector(N1MAX)
ccc FUNCTION
        REAL*8  Do3DScalarProduct
        REAL*8  GetDiscretizationResidual
ccc Passate per il calcolo del residual (temporal derivative)
        REAL*8  ET1E(NTEMAX),PEL(NTEMAX)
        REAL*8  PNEW(NMAX),PTIMEP(NMAX),SW_NODE(NMAX),SWP_NODE(NMAX)
        REAL*8  DELTAT
ccc passate per il l'assemblaggio dei flussi 
        INTEGER ISIDE(3,NFACEMAX)
        INTEGER SIDE_CNC(4,NTEMAX)
        REAL*8  AREFACE(NFACEMAX),Faces_ReferenceVector(3,NFACEMAX)
        REAL*8  Faces_FaceFlux(3,NFACEMAX)
        REAL*8  ElementCentroid(3,NTEMAX),FaceCentroid(3,NFACEMAX)

ccc Assemblo la parte relativa al flusso di Galerkin 
        DO J=1,Node_ElementCount(Node)
           CurrentElement=Node_ElementIDs(Node,J) 
ccc Faccio un ciclo sui nodi degli elementi 
           DO K=1,4
ccc Becco qual é in locale il nodo di riferimento
              IF(TETRA(K,CurrentElement).EQ.Node)THEN
                LocalNodeNumber=K
              END IF 
           END DO
ccccc Assemblo la parte relativa alla matrice di rigidezza 
           ResidualVector(J)=ResidualVector(J)+
     1     GetDiscretizationResidual(LocalNodeNumber,CurrentElement,
     1     TETRA,UU,VV,WW,X,Y,Z,VOLU)   
        END DO 
cc Adesso va assemblata la parte relativa alla discretizzazione 
cc temporale, sia di psi che di sw 
cc CARLOTTA CHECK IT LATER
        DO J=1,Node_ElementCount(Node)
           CurrentElement=Node_ElementIDs(Node,J)
           ResidualVector(J)=ResidualVector(J)-
     1     ET1E(CurrentElement)*0.25*VOLU(CurrentElement)*
     1     (PNEW(Node)-PTIMEP(Node))/DELTAT-
     1     PEL(CurrentElement)*0.25*VOLU(CurrentElement)*
     1     (SW_NODE(Node)-SWP_NODE(Node))/DELTAT
        END DO
ccc Sommo al residual la parte relativa ai flussi 
ccc che non conservono massa 
      DO J = 1,Node_ElementCount(Node)
        CurrentElement = Node_ElementIDs(Node,J)
        DO K = 1,4
              CurrentFaceID = SIDE_CNC(K,CurrentElement)
ccc Non voglio considerare la faccia che non é collegata al 
ccc nodo di riferimento. Se uno dei nodi della faccia é=a quello di 
ccc riferimento allora la faccia é da considerare
          IF ((ISIDE(1,CurrentFaceID).EQ.Node).OR.
     1        (ISIDE(2,CurrentFaceID).EQ.Node).OR.
     1        (ISIDE(3,CurrentFaceID).EQ.Node)) THEN
!           Get Current Face Area
            CurrentFaceArea = AREFACE(CurrentFaceID)
!           Get Face Normal
ccc Introduco in ciclo IF nel caso in cui la faccia é interna, allora 
ccc assemblo al residuals il flusso calcolato in precedenza come la
ccc media tra i flussi dei due elementi che condivide la faccia
ccc CARLOTTA CHECK
            CALL EvalElementFaceNormal(CurrentElement,K,
     1           Normal,Faces_ReferenceVector,XC,YC,ZC,
     1           SIDE_CNC,FaceCentroid)
c            IF(Faces_ElementCount(CurrentFaceID).EQ.2)THEN        
c            CALL EvalElementFaceNormal(CurrentElement,LoopC,
c     r                                 CurrentFaceNormal)
!           Make a local copy of the face flux
            DO L = 1,3
              LocalFaceFlux(L) =
     1        Faces_FaceFlux(L,CurrentFaceID)
            END DO
            FaceNormalFlux = Do3DScalarProduct(Normal,
     1                       LocalFaceFlux)
            ResidualVector(J) = ResidualVector(J)-(1.0d0/3.0d0)*
     1                          FaceNormalFlux*CurrentFaceArea
c           END IF
          END IF
        END DO
      END DO
c            ELSE
cccc Caso in cui la faccia é di Boundary 
c            IF(Node_Type(Node).EQ.0)THEN
cccc Inizializza i flussi sulle facce     
c              DO L=1,3
c               Faces_FaceFluxLocal(CurrentFace,L)=0.0d0
c              END DO     
c   Flux(1)=UU(CurrentElement)
c              Flux(2)=VV(CurrentElement)
c              FLux(3)=WW(CurrentElement)
c            Faces_FaceFluxLocal(CurrentFace,1)=UU(CurrentElement)
c            Faces_FaceFluxLocal(CurrentFace,2)=VV(CurrentElement)
c            Faces_FaceFluxLocal(CurrentFace,3)=WW(CurrentElement)
c            DO L=1,3
c             LocalFaceFlux(1) = UU(CurrentElement)
c             LocalFaceFlux(2) = VV(CurrentElement)
c             LocalFaceFlux(3) = WW(CurrentElement)
cc            END DO  
c             FaceNormalFlux = Do3DScalarProduct(Normal,
c     r                        LocalFaceFlux)
c             ResidualVector(J) = ResidualVector(J)-(1.0d0/3.0d0)*
c     r                           FaceNormalFlux*CurrentFaceArea
c            END IF
cc caso di Neumann Nullo 
c            IF(Node_Type(Node).EQ.2)THEN
cc attenta deve essere un flusso non volumetrico
c              CurrentNeumannFlux = Node_NeumannFlux(Node)
c                     CALL EvalElementFaceNormal(I,J,CurrentFaceNormal,
c     1                    Faces_ReferenceVector,ElementCentroid,
c     1                    Elements_FaceIDs,FaceCentroid)
c            DO L=1,3
c               LocalFaceFlux(L) = CurrentNeumannFlux*Normal(L)
c            END DO 
c            FaceNormalFlux = Do3DScalarProduct(Normal,
c     r                       LocalFaceFlux)
c            ResidualVector(J) = ResidualVector(J)-(1.0d0/3.0d0)*
c     r                          FaceNormalFlux*CurrentFaceArea

c            END IF
          
      RETURN
      END SUBROUTINE
ccc subroutine e function utlizzate per l'assemblaggio del residual
cccc per il Calcolo del Residuo relativo a un elemento proiettato 
cccc su un nodo 
      REAL*8 FUNCTION GetDiscretizationResidual(LocalNodeNumber,
     1                CurrentElement,TETRA,UU,VV,WW,X,Y,Z,VOLU)
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
!     External Functions
      REAL*8   Do3DScalarProduct
!     Arguments
      INTEGER LocalNodeNumber,CurrentElement
!     Local Variables
      REAL*8 ShapeDeriv(3)
      REAL*8 LocalFlux(3)
      REAL*8 IntProduct
ccc passate e usate qui
      INTEGER TETRA(5,NTEMAX) 
      REAl*8 UU(NTEMAX),VV(NTEMAX),WW(NTEMAX),VOLU(NTEMAX)
      REAL*8 X(NMAX),Y(NMAX),Z(NMAX)
       
!     Get Shape Function Derivative For the element
      CALL EvalTetraShapeDerivativeGlobal(CurrentElement,
     1                    LocalNodeNumber,TETRA,ShapeDeriv,X,Y,Z)
      LocalFlux(1)=UU(CurrentElement)
      LocalFlux(2)=VV(CurrentElement)
      LocalFlux(3)=WW(CurrentElement)
      IntProduct = Do3DScalarProduct(LocalFlux,ShapeDeriv)
      GetDiscretizationResidual = IntProduct*
     1       VOLU(CurrentElement)
      RETURN
      END FUNCTION

cccc  Calcolare Grad delle funzioni test 
      SUBROUTINE EvalTetraShapeDerivativeGlobal(CurrentElement,
     1           LocalNodeNumber,TETRA,ShapeDeriv,X,Y,Z)
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
!     Arguments
      INTEGER CurrentElement,LocalNodeNumber
cc passate da passare 
      INTEGER TETRA(5,NTEMAX)
      REAL*8 ShapeDeriv(3)
      REAL*8 X(NMAX),Y(NMAX),Z(NMAX)
!     Local Variables
      INTEGER I,J
      REAL*8 LocalDerivatives(3)
      REAL*8 JacobianMatrix(3,3)
!     Eval Local Shape Derivatives
      CALL EvalTetraShapeDerivativeLocal(LocalNodeNumber,
     r                                    LocalDerivatives)
!     Eval Element Jacobian Matrix
      CALL EvalTetraJacobianMatrix(CurrentElement,JacobianMatrix,TETRA,
     1     X,Y,Z)
!     Multiply
      DO I = 1,3
        ShapeDeriv(I) = 0.0d0
        DO J = 1,3
          ShapeDeriv(I) = ShapeDeriv(I)+
     r    JacobianMatrix(I,J)*LocalDerivatives(J)
        END DO
      END DO
      END SUBROUTINE
ccccc
      SUBROUTINE EvalTetraShapeDerivativeLocal(LocalNodeNumber,
     1                                          ShapeLocalDerivatives)
      IMPLICIT NONE
!     Arguments
      INTEGER LocalNodeNumber
      REAL*8 ShapeLocalDerivatives(3)
      IF (LocalNodeNumber.EQ.1) THEN
!        d_N1/d_Psi
         ShapeLocalDerivatives(1) = 1.0d0
!        d_N1/d_Eta
         ShapeLocalDerivatives(2) = 0.0d0
!        d_N1/d_Chi
         ShapeLocalDerivatives(3) = 0.0d0
      ELSE IF (LocalNodeNumber.EQ.2) THEN
!        d_N2/d_Psi
         ShapeLocalDerivatives(1) = 0.0d0
!        d_N2/d_Eta
         ShapeLocalDerivatives(2) = 1.0d0
!        d_N2/d_Chi
         ShapeLocalDerivatives(3) = 0.0d0
      ELSE IF (LocalNodeNumber.EQ.3) THEN
!        d_N3/d_Psi
         ShapeLocalDerivatives(1) = 0.0d0
!        d_N3/d_Eta
         ShapeLocalDerivatives(2) = 0.0d0
!        d_N3/d_Chi
         ShapeLocalDerivatives(3) = 1.0d0
      ELSE IF (LocalNodeNumber.EQ.4) THEN
!        d_N4/d_Psi
         ShapeLocalDerivatives(1) = -1.0d0
!        d_N4/d_Eta
         ShapeLocalDerivatives(2) = -1.0d0
!        d_N4/d_Chi
         ShapeLocalDerivatives(3) = -1.0d0
      END IF
      END SUBROUTINE
cccccc
      SUBROUTINE EvalTetraJacobianMatrix(CurrentElement,JMat,TETRA,
     1           X,Y,Z)
      IMPLICIT NONE
!     Include Global Variables
      INCLUDE 'CATHY.H'
!     Arguments
      INTEGER CurrentElement,TETRA(5,NTEMAX)
      REAL*8 JMat(3,3)
!     Local Variables
      INTEGER I
      REAL*8 dX_dPsi,dX_dEta,dX_dChi
      REAL*8 dY_dPsi,dY_dEta,dY_dChi
      REAL*8 dZ_dPsi,dZ_dEta,dZ_dChi
      INTEGER CurrentNode
      REAL*8 CurrentNodeDerivs(3)
      REAL*8 Mat(3,3)
cccc Passate 
      REAL*8 X(NMAX),Y(NMAX),Z(NMAX)
!     Initialize
      dX_dPsi = 0.0d0
      dX_dEta = 0.0d0
      dX_dChi = 0.0d0
      dY_dPsi = 0.0d0
      dY_dEta = 0.0d0
      dY_dChi = 0.0d0
      dZ_dPsi = 0.0d0
      dZ_dEta = 0.0d0
      dZ_dChi = 0.0d0
      DO I = 1,4
!       Eval Current Node Derivs
        CurrentNode = TETRA(I,CurrentElement)
        CALL EvalTetraShapeDerivativeLocal(I,CurrentNodeDerivs)
        dX_dPsi = dX_dPsi+CurrentNodeDerivs(1)*
     r            X(CurrentNode)
        dX_dEta = dX_dEta+CurrentNodeDerivs(2)*
     r            X(CurrentNode)
        dX_dChi = dX_dChi+CurrentNodeDerivs(3)*
     r            X(CurrentNode)
        dY_dPsi = dY_dPsi+CurrentNodeDerivs(1)*
     r            Y(CurrentNode)
        dY_dEta = dY_dEta+CurrentNodeDerivs(2)*
     r            Y(CurrentNode)
        dY_dChi = dY_dChi+CurrentNodeDerivs(3)*
     r            Y(CurrentNode)
        dZ_dPsi = dZ_dPsi+CurrentNodeDerivs(1)*
     r            Z(CurrentNode)
        dZ_dEta = dZ_dEta+CurrentNodeDerivs(2)*
     r            Z(CurrentNode)
        dZ_dChi = dZ_dChi+CurrentNodeDerivs(3)*
     r            Z(CurrentNode)
      END DO
!     Assign Matrix Coeffs
!     1
      Mat(1,1) = dX_dPsi
      Mat(1,2) = dY_dPsi
      Mat(1,3) = dZ_dPsi
!     2
      Mat(2,1) = dX_dEta
      Mat(2,2) = dY_dEta
      Mat(2,3) = dZ_dEta
!     3
      Mat(3,1) = dX_dChi
      Mat(3,2) = dY_dChi
      Mat(3,3) = dZ_dChi
!     Allocate
      CALL Invert3x3Matrix(Mat,JMat)
      END SUBROUTINE

cccc Invertire matrice 3x3
      SUBROUTINE Invert3x3Matrix(Mat,Inverse)
      IMPLICIT NONE
!     Arguments
      REAL*8 Mat(3,3),Inverse(3,3)
!     Local Variables
      REAL*8 DetJ
!     Valuta il Determinante
      DetJ = Mat(1,1)*(Mat(2,2)*Mat(3,3)-Mat(2,3)*Mat(3,2))-
     r       Mat(1,2)*(Mat(2,1)*Mat(3,3)-Mat(2,3)*Mat(3,1))+
     r       Mat(1,3)*(Mat(2,1)*Mat(3,2)-Mat(2,2)*Mat(3,1))
!     Write Matrix
      Inverse(1,1) = (1.0d0/DetJ)*(Mat(2,2)*Mat(3,3)-Mat(2,3)*Mat(3,2))
      Inverse(2,2) = (1.0d0/DetJ)*(Mat(1,1)*Mat(3,3)-Mat(1,3)*Mat(3,1))
      Inverse(3,3) = (1.0d0/DetJ)*(Mat(1,1)*Mat(2,2)-Mat(1,2)*Mat(2,1))
      Inverse(1,2) =-(1.0d0/DetJ)*(Mat(1,2)*Mat(3,3)-Mat(1,3)*Mat(3,2))
      Inverse(1,3) = (1.0d0/DetJ)*(Mat(1,2)*Mat(2,3)-Mat(2,2)*Mat(1,3))
      Inverse(2,1) =-(1.0d0/DetJ)*(Mat(2,1)*Mat(3,3)-Mat(2,3)*Mat(3,1))
      Inverse(2,3) =-(1.0d0/DetJ)*(Mat(1,1)*Mat(2,3)-Mat(1,3)*Mat(2,1))
      Inverse(3,1) = (1.0d0/DetJ)*(Mat(2,1)*Mat(3,2)-Mat(3,1)*Mat(2,2))
      Inverse(3,2) =-(1.0d0/DetJ)*(Mat(1,1)*Mat(3,2)-Mat(1,2)*Mat(3,1))
      END SUBROUTINE
c-----------FINE parte relativa al residual--------------------------- 
ccc---------CALCOLO DELLA MATRICE DEI COEFFICIENTI LOCALE------------- 
ccc Subroutine principale per il calcolo della matrice
      SUBROUTINE AssembleMatrix(Node,LocMat,AREFACE,XC,YC,ZC,
     1           Node_ElementIDs,Node_ElementCount,
     1           SIDE_CNC,Node_FaceIDs,Node_FaceCount,
     1           PLIST,Faces_FaceType)

      IMPLICIT NONE 
      INCLUDE 'CATHY.H'
cc Locali
      INTEGER J,K,L
      INTEGER CurrentElementID,LocalNodeFace(MaxFConToNode)
      INTEGER CurrentFaceID
      INTEGER LocalNodeElements(N1MAX)
      INTEGER OtherElementID,OtherElementIDx
      LOGICAL IsInNodeFaceList
cc Passate 
      INTEGER Node
      INTEGER Node_ElementIDs(NMAX,N1MAX),Node_ElementCount(NMAX)
      INTEGER Node_FaceIDs(NMAX,MaxFConToNode),Node_FaceCount(NMAX)
      INTEGER SIDE_CNC(4,NTEMAX)
      INTEGER PLIST(2,NFACEMAX),Faces_FaceType(NFACEMAX)
      REAL*8  LocMat(N1MAX,N1MAX)
      REAl*8  AREFACE(NFACEMAX)
      REAl*8  XC(NTEMAX),YC(NTEMAX),ZC(NTEMAX)

      DO J = 1,Node_ElementCount(Node)
         CurrentElementID = Node_ElementIDs(Node,J)
!        Eval Current Element Centroid
!        Edge Loop
         DO K = 1,4
           CurrentFaceID =
     1     SIDE_CNC(K,CurrentElementID)
!          Valuta la Lunghezza dello spigolo
c          LocalFaceArea = EvalFaceArea(CurrentFaceID)
!          Fai una Copia locale della lista spigoli
           DO L = 1,Node_FaceCount(Node)
              LocalNodeFace(L) = Node_FaceIDs(Node,L)
           END DO
!          Se Lo spigolo fa parte della stella
           IF (IsInNodeFaceList(CurrentFaceID,
     1                          Node_FaceCount(Node),
     1                          LocalNodeFace)) THEN
!            Trova l'altro elemento
             IF (PLIST(1,CurrentFaceID).EQ.
     1           CurrentElementID) THEN
                 OtherElementID = PLIST(2,CurrentFaceID)
             ELSE
              OtherElementID = PLIST(1,CurrentFaceID)
             END IF
!           Se è connesso qualche altro elemento
            IF (OtherElementID.GT.0) THEN
!             Fai una copia locale della lista di elementi
              DO L = 1,Node_ElementCount(Node)
                LocalNodeElements(L) = Node_ElementIDs(Node,L)
              END DO
!             Trova l'indice locale dell'elemento
cc
              CALL GetFromNodeElements(OtherElementID,
     1                   Node_ElementCount(Node),
     1                   LocalNodeElements,
     1                   OtherElementIdx)
              LocMat(J,J) = LocMat(J,J)+
     1                  (1.0d0/3.0d0)*AREFACE(CurrentFaceID)
              LocMat(J,OtherElementIdx) = 
     1        LocMat(J,OtherElementIdx)
     1        -(1.0d0/3.0d0)*AREFACE(CurrentFaceID)
            ELSE
cccc caso in cui faccia di Dirichlet
cccc Se il nodo é di Dirichlet, la faccia é di boundary, allora la
cccc faccia é di Dirichlet
              IF(Faces_FaceType(CurrentFaceID).EQ.0)THEN
c              IF (Faces_FaceType(CurrentFaceID).EQ.
c     1            0) THEN
!               Assembla Solo il Diagonale
                LocMat(J,J) = LocMat(J,J) +
     1          (1.0d0/3.0d0)*AREFACE(CurrentFaceID)
              END IF
            END IF 
         END IF
        END DO 
      END DO 
      RETURN 

      END SUBROUTINE
cccc Logical function per valutare quale faccia é quella considerata 
      LOGICAL FUNCTION IsInNodeFaceList(ID,Count,List)
      IMPLICIT NONE
!     Include Global Variables
      INCLUDE 'CATHY.H'
!     Variabili Locali
      INTEGER ID,I
      INTEGER Count
      INTEGER List(MaxFConToNode)
!     Corpo Procedura
      IsInNodeFaceList = .FALSE.
      DO I = 1,Count
        IF (List(I).EQ.ID) THEN
          IsInNodeFaceList = .TRUE.
          RETURN
        END IF
      END DO
      END FUNCTION

cccccccc
! ---------------------------------------------------------------------
! RECUPERA ELEMENTI AFFERENTI AD UN DATO NODO  CARLOTTA OKOKOKO
! ---------------------------------------------------------------------
      SUBROUTINE GetFromNodeElements(ID,Count,List,Position)
      IMPLICIT NONE
!     Include Global Variables
      INCLUDE 'CATHY.H'
!     Variabili Locali
      INTEGER ID,I
      INTEGER Count
      INTEGER List(N1MAX)
      INTEGER Position
      DO I = 1,Count
        IF (List(I).EQ.ID) THEN
          Position = I
          RETURN
        END IF
      END DO
      END SUBROUTINE
c-------------------- FINE SUBROUTINE per IL CALCOLO 
c------------- Risolvere il sistema locale ----------------c
c-----------SOLUTION OF THE SYSTEM-------------------------------------
      SUBROUTINE SolveFor3DAuxiliaryFluxes(CurrentNode,
     1           LocMat,ResidualVector,Node_ElementCount,
     1           Node_FaceCount,PLIST,Node_FaceIDs)
      IMPLICIT NONE
!     Include Global Variables
      INCLUDE 'CATHY.H'
!     External Functions
      LOGICAL IsInternal3DMeshNode
      INTEGER GaussElimination,I
!      LOGICAL IsNodeInFace
!     Arguments
c      LOGICAL FUNCTION IsInternal3DMeshNode(CurrentNode,Node_FaceCount,
c     1                 Faces_ElementCount,Node_FaceIDs)
      INTEGER Node_ElementCount(NMAX)
      INTEGER CurrentNode,CurrentFaceID,cont
      INTEGER Node_FaceCount(NMAX),PLIST(2,NFACEMAX)
      INTEGER Node_FaceIDs(NMAX,MaxFConToNode)
      REAL*8  LocMat(N1MAX,N1MAX)
      REAL*8  ResidualVector(N1MAX)
      REAL*8  MathZero
!     Local Variables
!      INTEGER*4 LoopA,LoopB
!      INTEGER*4 CurrentElement,CurrentFaceID
      INTEGER GaussError
!      LOGICAL FoundNeumann
      REAL*8  Residual
      REAL*8  New_Mat(N1MAX,N1MAX)
      REAL*8  New_RHS(N1MAX)
      INTEGER Count
      MathZero=1.0E-07
!     Function Core
      Residual = -1.0d0
c         cont=0
c      DO I=1,Node_FaceCount(CurrentNode)
c         CurrentFaceID=Node_FaceIDs(CurrentNode,I)
c         IF(Faces_ElementCount(CurrentFaceID).EQ.2)THEN
c           cont=cont+1
c         END IF 
c      END DO 
c      IF(CurrentNode.Eq.1)THEN
c      write(*,*) cont,Node_FaceCount(1)
c      DO I=1,Node_FaceCount(1)
c      write(*,*)Node_FaceIDs(1,I),
c     1     Faces_ElementCount(Node_FaceIDs(1,I))
c      END DO 
c      END IF
c      STOP
c      write(*,*) CurrentNode,Node_faceCount(CurrentNode),cont
c      write(*,*) 'main1'
c       IF(cont.eq.Node_FaceCount(CurrentNode))THEN
      IF (IsInternal3DMeshNode(CurrentNode,Node_FaceCount,
     1    PLIST,Node_FaceIDs)) THEN
!       Nodi Interni, Elimina la costante
        CALL ModifyLocalMat(LocMat,ResidualVector,
     1                      Node_ElementCount(CurrentNode),1,0.0d0)
!       Solve Equation Set
        GaussError = 
     1  GaussElimination(Node_ElementCount(CurrentNode),
     1                   LocMat,ResidualVector,Residual,
     1                   CurrentNode)
        IF (GaussError.NE.0) THEN
          WRITE(*,*)'Error During Gauss Elimination for Internal Node.'
          STOP
        END IF
      ELSE
        IF ((Node_ElementCount(CurrentNode).EQ.1).AND.
     r     (ABS(LocMat(1,1)).LT.MathZero)) THEN
            ResidualVector(1) = 0.0d0
        ELSE
!         Elimina Tutte le equazioni con spigoli alla Neumann
!          DO LoopA = 1,MeshNodes_ElementCount(CurrentNode)
!            CurrentElement = MeshNodes_ElementIDs(CurrentNode,LoopA)
!           Controlla se vo sono Condizioni di Neumann Applicate
!            FoundNeumann = .FALSE.
!            DO LoopB = 1,MeshElements_EdgeFaceCount(CurrentEle/ment)
!              CurrentFaceID =
!     r        MeshElements_EdgeFaceIDs(CurrentElement,LoopB)
!     r        .AND.(IsNodeInFace(CurrentNode,CurrentFaceID))) THEN
!                FoundNeumann = .TRUE.
!              END IF
!            END DO
!           If Found Set the Correction equal to 0.0
!            IF (FoundNeumann) THEN
!              CALL ModifyLocalMat(Local_Mat,RHS_Vector,
!     r             MeshNodes_ElementCount(CurrentNode),LoopA,0.0d0)
!            END IF
!          END DO
!         Copy
          CALL CopyMatrix(Node_ElementCount(CurrentNode),
     r                    LocMat,New_Mat)
          CALL CopyVector(Node_ElementCount(CurrentNode),
     r                    ResidualVector,New_RHS)
          GaussError =
     r    GaussElimination(Node_ElementCount(CurrentNode),New_Mat,
     r                     New_RHS,Residual,CurrentNode)
          Count = 0
          DO WHILE (GaussError.NE.0)
            Count = Count + 1
           CALL ModifyLocalMat(LocMat,ResidualVector,
     r           Node_ElementCount(CurrentNode),Count,0.0d0)
            CALL CopyMatrix(Node_ElementCount(CurrentNode),
     r                      LocMat,New_Mat)
           CALL CopyVector(Node_ElementCount(CurrentNode),
     r                      ResidualVector,New_RHS)
!           Solve Equation Set
            GaussError =
     r      GaussElimination(Node_ElementCount(CurrentNode),
     r                       New_Mat,New_RHS,Residual,CurrentNode)
            END DO
!         Copy Solution
          CALL CopyVector(Node_ElementCount(CurrentNode),
     r                    New_RHS,ResidualVector)
        END IF
      END IF
      END SUBROUTINE
! ---------------------------------------------------------------------
      LOGICAL FUNCTION IsInternal3DMeshNode(CurrentNode,Node_FaceCount,
     1                 PLIST,Node_FaceIDs)
      IMPLICIT NONE
!     Include Global Variables
      INCLUDE 'CATHY.H'
!     Arguments
      INTEGER CurrentNode
      INTEGER Node_FaceCount(NMAX),PLIST(2,NFACEMAX)
      INTEGER Node_FaceIDs(NMAX,MaxFConToNode)
!     Local Variables
      INTEGER I
      LOGICAL LocalCheck
      INTEGER CurrentFaceID
!     Function Core
      LocalCheck = .TRUE.
      DO I = 1,Node_FaceCount(CurrentNode)
        CurrentFaceID = Node_FaceIDs(CurrentNode,I)
        LocalCheck = (LocalCheck.AND.
     r               (PLIST(2,CurrentFaceID).NE.0))
      END DO

      IsInternal3DMeshNode = LocalCheck
      RETURN
      END FUNCTION
c-----------------------------------------------------------------------------
c----------------------Modifica della matrice locale--------------------------
c-----------------------------------------------------------------------------      
      SUBROUTINE ModifyLocalMat(LocMat,ResidualVector,
     r                          Count,Index,RHSValue)
      IMPLICIT NONE
!     Include Global Variables
      INCLUDE 'CATHY.H'
!     Argomenti
      INTEGER Count
      INTEGER Index
      REAL*8  LocMat(N1MAX,N1MAX)
      REAL*8  ResidualVector(N1MAX)
      REAL*8  RHSValue
!     Variabili Locali
      INTEGER I
!     Set the RHS equal to zero
      ResidualVector(Index) = RHSValue
!     Set the row with the diagonal equal to one and zero elsewere
      DO I = 1,Count
        IF (I.EQ.Index) THEN
          LocMat(Index,I) = 1.0d0
        ELSE
          LocMat(Index,I) = 0.0d0
        END IF
      END DO
      END SUBROUTINE
! RISOLVI UN SISTEMA LINEARE CON IL METODO DI GAUSS CARLOTTA OKOKOK
! ---------------------------------------------------------------------
      INTEGER FUNCTION GaussElimination(Size,AMat,Rhs,ResidualNorm,
     1               CurrentNode)
      IMPLICIT NONE
!     Include Global Variables
      INCLUDE 'CATHY.H'
!     Argomenti
      INTEGER Size,CurrentNode
      REAl*8 MathZero
      REAL*8 AMat(N1MAX,N1MAX)
      REAL*8 Rhs(N1MAX)
      REAL*8 ResidualNorm
!     Local Variables
      INTEGER A,B,C
      REAL*8 InitialNorm
      REAL*8 Factor
      REAL*8 Solution(Size)
      REAL*8 Residual(Size)
!     Functions and Subs
      REAL*8 EucNorm
      MathZero=1.0E-07
      GaussElimination = 0
      InitialNorm = EucNorm(Size,Rhs);
!     Reduction To Triangular Form
      DO A = 1,Size
        DO B = (A+1),Size
          IF (ABS(AMat(A,A)).LT.MathZero) THEN
            GaussElimination = 1
            RETURN
          ELSE
            Factor = -(AMat(B,A)/AMat(A,A));
          END IF
          DO C = 1,Size
            AMat(B,C) = AMat(B,C) + Factor*AMat(A,C)
          END DO
!     RHS
          Rhs(B) = Rhs(B) + Factor*Rhs(A)
        END DO
      END DO
!     Back Sub
      A = Size+1
      DO WHILE (A>1)
        A = A-1
        Solution(A) = Rhs(A)
        DO B = A+1,Size
          Solution(A) = Solution(A)-AMat(A,B)*Solution(B)
        END DO
        IF (ABS(AMat(A,A)).LT.MathZero) THEN
          GaussElimination = 1
          RETURN
        ELSE
          Solution(A) = Solution(A)/AMat(A,A)
        END IF
      END DO
!     Verify Error
      DO A = 1,Size
        Residual(A) = Rhs(A)
        DO B = 1,Size
          Residual(A) = Residual(A)-AMat(A,B)*Solution(B)
        END DO
      END DO
!     Copy the Solution Vector
      DO A = 1,Size
        Rhs(A) = Solution(A)
      END DO
      IF (ABS(InitialNorm).LT.MathZero) THEN
        ResidualNorm = 0.0d0
      ELSE
        ResidualNorm = (EucNorm(Size,Residual)/InitialNorm);
      END IF
      RETURN
      END FUNCTION
cc---- Norma Euclidea-------------------------------
      REAL*8 FUNCTION EucNorm(Size,Rhs)
      IMPLICIT NONE
!     Include Global Variables
      INCLUDE 'CATHY.H'
!     Arguments
      INTEGER Size
      REAL*8 Rhs(N1MAX)
!     Local Variales
      INTEGER I
      EucNorm = 0.0d0
      DO I = 1,Size
         EucNorm = EucNorm + Rhs(I) * Rhs(I)
      END DO
      EucNorm = SQRT(EucNorm)
      END FUNCTION
c----------------------------------------------------------------------
c-----------------------------COPYMATRIX---------------------------------
c-------------------------------------------------------------------------
      SUBROUTINE CopyMatrix(Size,OrigMat,TargetMat)
      IMPLICIT NONE
!     Include Global Variables
      INCLUDE  'CATHY.H'
      INTEGER  I,J
      INTEGER  Size
      REAL*8   OrigMat(N1MAX,N1MAX)
      REAL*8 TargetMat(N1MAX,N1MAX)
      DO I = 1,Size
        DO J = 1,Size
          TargetMat(I,J) = OrigMat(I,J)
        END DO
      END DO
      END SUBROUTINE
!-----------------------------------------------------------------------
      SUBROUTINE CopyVector(Size,OrigVector,TargetVector)
      IMPLICIT NONE
!     Include Global Variables
      INCLUDE 'CATHY.H'
      INTEGER  I
      INTEGER  Size
      REAL*8 OrigVector(N1MAX)
      REAL*8 TargetVector(N1MAX)
      DO I = 1,Size
         TargetVector(I) = OrigVector(I)
      END DO
      END SUBROUTINE
!-----------------------------------------------------------------------
