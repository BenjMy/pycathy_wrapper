       SUBROUTINE VELOCITY_MAIN(N,X,Y,Z,PNEW,PTIMEP,SW_NODE,SWP_NODE,
     1                       DELTAT,  
     1                       NT,TETRA,VOLU,IVOL,ET1E,PEL,UU,VV,WW,
     1                       Node_ElementCount,Node_ElementIDs,ISIDE,
     1                       PLIST,FaceCentroid,
     1                       AREFACE,Faces_NeumannFlux,
     1                       Faces_FaceFlux,Node_FaceCount,
     1                       Node_FaceIDs,Faces_FaceType,NFACE,
     1                       Faces_ReferenceVector,SIDE_CNC,
     1                       XC,YC,ZC,time)

       IMPLICIT NONE 

       INCLUDE 'CATHY.H'
       INTEGER CurrentFace
cc variabili locali
       INTEGER I,J
       INTEGER cont,Node1,Node2,Node3,CurrentElement,LocalFace
ccc CARLOTTA TEMPORANEE
       INTEGER Element1,Element2
       REAL*8  Normal1(3),Normal2(3)
       REAL*8  V1,V2
       REAL*8  Error(NFACEMAX)
       REAL*8  FaceFlux(3)
       REAL*8  Do3DScalarProduct
       INTEGER Faccia
cc INPUT 
       REAL*8  time
cc INPUT Per nodi 
       INTEGER N
       REAL*8  X(NMAX),Y(NMAX),Z(NMAX)
       REAL*8  PNEW(NMAX),PTIMEP(NMAX),SW_NODE(NMAX),SWP_NODE(NMAX)
cc Per elementi 
       INTEGER NT,TETRA(5,NTEMAX)
       INTEGER IVOL(NTEMAX)
       REAL*8  VOLU(NTEMAX),ET1E(NTEMAX),PEL(NTEMAX)
       REAL*8  UU(NTEMAX),VV(NTEMAX),WW(NTEMAX)
cc per il calcolo del residual in derivata temporale 
       REAl*8  DELTAT
ccc Variabili introdotte per velocity reconstruction 
ccc NODI 
       INTEGER Node_ElementCount(NMAX),Node_ElementIDs(NMAX,N1MAX)
       INTEGER Faces_FaceType(NFACEMAX)
ccc FACCE
       INTEGER NFACE
       INTEGER ISIDE(3,NFACEMAX),PLIST(2,NFACEMAX)
       INTEGER Node_FaceIDs(NMAX,MaxFConToNode)
       INTEGER Node_FaceCount(NMAX)
       REAL*8  Faces_ReferenceVector(3,NFACEMAX)
       REAL*8  FaceCentroid(3,NFACEMAX),AREFACE(NFACEMAX)
ccc Generati qui dentro 
       REAL*8  Faces_FaceFlux(3,NFACEMAX),Faces_NeumannFlux(NFACEMAX)
ccc ELEMENTI
       INTEGER SIDE_CNC(4,NTEMAX)
       REAL*8  XC(NTEMAX),YC(NTEMAX),ZC(NTEMAX)
ccc FUNCTION Per il calcolo dell'IMBALANCE 
       REAL*8  EvalNetElementFlux
       REAL*8  CurrentImbalancePrima,CurrentImbalanceDopo
cccc SOLUZIONE che arriva da applyLN
       REAl*8 Faces_FluxCorrection(NFACEMAX,3)
       REAl*8 Elements_FluxXCentroid(NTEMAX)
       REAl*8 Elements_FluxYCentroid(NTEMAX)
       REAL*8 Elements_FluxZCentroid(NTEMAX)
          
       CALL VELOCITYINITAL(NFACE,NT,UU,VV,WW,
     1      SIDE_CNC,PLIST,
     1      Faces_FaceType, 
     1      ISIDE,
     1      XC,YC,ZC,FaceCentroid,
     1      Faces_ReferenceVector,
     1      Faces_NeumannFlux,
     1      Faces_FaceFlux)
c       DO I=1,Node_FaceCount(63)
c        CurrentFace=Node_FaceIDs(63,I)
c          write(202,*) CurrentFace, Faces_NodeIDs(CurrentFace,1),
c     1    Faces_NodeIDs(CurrentFace,2),Faces_NodeIDs(CurrentFace,3)
c       END DO 

ccccc CALCOLO degli squilibri negli elementi prima del post-processing
ccccc come somma del current imbalance di ciascun elemento
cc Imbalance calculated by considering element by element the internal
cc faces 
       CurrentImbalancePrima=0.0d0
      DO I=1,NT
        CurrentImbalancePrima = CurrentImbalancePrima + 
     1          EvalNetElementFlux(I,PLIST, 
     1          SIDE_CNC,AREFACE,Faces_ReferenceVector,
     1          XC,YC,ZC,FaceCentroid,Faces_FaceFlux,
     1          UU,VV,WW)
       END DO
CCC CARLOTTA Calcolo un altro imbalance qui
CC DA CANCELLARE
c       IF(TIME.GT.2000.)THEN
c       DO I=1,NFACE
c         IF(PLIST(2,I).NE.0)THEN
cc  Calcolo velocità secondo elemento PLIST(1,FACCIA)
c          Element1=PLIST(1,I)
c          Element2=PLIST(2,I)
c          CALL EvalElementFaceNormal2(Element1,I,Normal1,
c     1         Faces_ReferenceVector,XC,YC,ZC,FaceCentroid)
c          CALL EvalElementFaceNormal2(Element2,I,Normal2,
c     1         Faces_ReferenceVector,XC,YC,ZC,FaceCentroid)
c          V1=UU(Element1)*Normal1(1)+VV(Element1)*Normal1(2)+
c     1    WW(Element1)*Normal1(3)
c          V2=UU(Element2)*Normal2(1)+
c     1    VV(Element2)*Normal2(2)+WW(Element2)*Normal2(3)
c          Error(I)=V1+V2
c          Write(1110,*)I,V1,V2,Error(I)
c         ELSE
c          write(1110,*)I,'Faccia di bordo'
c         END IF
c       END DO
c       END IF 

ccc Applico al correzione di Larson Niklasson. Tutte le subroutine 
ccc Relative al calcolo le trovo in larson-niklasson.f
       CALL ApplyLN(N,NT,Node_ElementCount,Node_ElementIDs,
     1      TETRA,UU,VV,WW,X,Y,Z,VOLU,ET1E,PEL,PNEW,
     1      PTIMEP,SW_NODE,SWP_NODE,DELTAT,
     1      ISIDE,AREFACE,
     1      Faces_ReferenceVector,Faces_FaceFlux,
     1      XC,YC,ZC,FaceCentroid,
     1      SIDE_CNC,
     1      Faces_FaceType,Node_FaceIDs,
     1      Node_FaceCount,PLIST,
     1      Faces_FluxCorrection,
     1      Faces_NeumannFlux)
cccc Aggiorna i flussi 
       CALL UpdateFaceFluxes(NFACE,Faces_FaceFlux,
     1      Faces_FluxCorrection,Faces_ReferenceVector,
     1      PLIST,
     1      N,Node_FaceCount,Node_FaceIDs,
     1      SIDE_CNC,
     1      XC,YC,ZC,FaceCentroid,UU,VV,WW)
       CALL EvalCentroidElementFluxes(NT, 
     1           AREFACE,XC,YC,ZC,VOLU,TETRA,
     1           SIDE_CNC,
     1           Elements_FluxXCentroid,Elements_FluxYCentroid,
     1           Elements_FluxZCentroid,Faces_ReferenceVector,
     1           FaceCentroid,X,Y,Z,ISIDE,Faces_FaceFlux)
      
        CurrentImbalanceDopo=0.0d0 
        DO I = 1,NT
        CurrentImbalanceDopo = CurrentImbalanceDopo+
     1          EvalNetElementFlux(I,PLIST,
     1          SIDE_CNC,AREFACE,Faces_ReferenceVector,
     1          XC,YC,ZC,FaceCentroid,Faces_FaceFlux,
     1          Elements_FluxXCentroid,Elements_FluxYCentroid, 
     1          Elements_FLuxZCentroid)
c          WRITE(98,1001)I,CurrentImbalance
c          IF (ABS(CurrentImbalance).GT.ABS(MaxImbalanceAfterLN)) THEN
c            MaxImbalanceAfterLN = CurrentImbalance
c          END IF
        END DO
c       cont=0 
c       CALL VTKRIS3D(NT,N,500+cont,TETRA,TIME,
c     1      PNEW,SWNEW,Elements_FluxXCentroid,Elements_FluxYCentroid,
c     1      Elements_FluxZCentroid,X,Y,Z,KS)
c       cont=cont+1
CCC CARLOTTA Calcolo un altro imbalance qui
c       IF(Time.GT.2000.)THEN
c       DO I=1,NFACE
c         IF(PLIST(2,I).NE.0)THEN
cc  Calcolo velocità secondo elemento PLIST(1,FACCIA)
c          Element1=PLIST(1,I)
c          Element2=PLIST(2,I)
c          CALL EvalElementFaceNormal2(Element1,I,Normal1,
c     1         Faces_ReferenceVector,XC,YC,ZC,FaceCentroid)
c          CALL EvalElementFaceNormal2(Element2,I,Normal2,
c     1         Faces_ReferenceVector,XC,YC,ZC,FaceCentroid)
c          V1=Elements_FluxXCentroid(Element1)*Normal1(1)+
c     1    Elements_FluxYCentroid(Element1)*Normal1(2)+
c     1    Elements_FluxZCentroid(Element1)*Normal1(3)
c          V2=Elements_FluxXCentroid(Element2)*Normal2(1)+
c     1    Elements_FLuxYCentroid(Element2)*Normal2(2)+
c     1    Elements_FluxZCentroid(Element2)*Normal2(3)
c          Error(I)=V1+V2
c          Write(1111,*)I,V1,V2,Error(I)
c         ELSE
c          write(1111,*)I,'Faccia di bordo'
c         END IF
c       END DO
c       STOP
c       END IF 

c      write(300,*)
c    1    Time,'PRIMA',CurrentImbalancePrima,'DOPO',CurrentImbalanceDopo
c       Print 2D Results
c     questa devi completarla
c        CALL Print3DResults
c     Scrivi un messaggio finale
c      WRITE(*,*) 'Esecuzione Completata!'
c      WRITE(*,1000)MaxImbalanceBeforeLN,MaxImbalanceAfterLN
c     Format Statements
c1000  FORMAT( 'Imbalance Before LN:',2X,E10.3,/,
c     r        'Imbalance After LN: ',2X,E10.3)
c1001  FORMAT(I8,E15.5)
c1002  FORMAT(I8,3E15.5)
c       WRITE(IOUT60) TIME 
c       WRITE(IOUT60) (Elements_FluxXCentroid(I),I=1,NT)
c       WRITE(IOUT60) (Elements_FluxYCentroid(I),I=1,NT)
c       WRITE(IOUT60) (Elements_FluxZCentroid(I),I=1,NT)
       DO I=1,NT
          UU(I)=Elements_FluxXCentroid(I)
          VV(I)=Elements_FluxYCentroid(I)
          WW(I)=Elements_FluxZCentroid(I)
       END DO  


       RETURN
       END 
