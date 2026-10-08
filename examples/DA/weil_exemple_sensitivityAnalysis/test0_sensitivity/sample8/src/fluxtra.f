C
C**************************  FLUXTRA ***********************************
C
C  compute Darcy velocity values for first time step
c  velocity is VN (normal) component
C  3D 
C***********************************************************************
C
      SUBROUTINE FLUXTRA(N,NFACE,NT,VELREC,
     1           PLIST,ISIDE,SIDE_CNC,
     2           Node_FaceCount,Node_FaceIDS,          
     3           X,Y,Z,XC,YC,ZC,Faces_ReferenceVector,FaceCentroid,
     4           AREFACE,VX,VY,VZ,FACEFLUX,
     5           Faces_FaceType,Faces_NeumannFlux,
     6           BKCFLOWINT_NODE,ARENOD,VN)
C
      IMPLICIT  NONE
      INCLUDE  'CATHY.H'
C     LOCAL
      INTEGER  I,K,j,itetra,jtetra,ltetra,iel1,iel2,iface
      INTEGER  Node1,Node2,Node3,Faccia
      REAL*8   vxface,vyface,vzface
      REAL*8   alpha, beta, gamma, prodscal
      REAL*8   xl(3), yl(3), zl(3)
      REAL*8   Normal1(3)
      REAL*8   FLUX(3),AREA(NMAX)
      REAL*8   Do3DScalarProduct
C     PASSED
      INTEGER  n,nface,nt,velrec
      INTEGER  plist(2,*),side_cnc(4,*)
      INTEGER  iside(3,*)  
      INTEGER  Faces_FaceType(*) 
      INTEGER  Node_faceCount(*),Node_FaceIDS(NMAX,*)
      real*8   x(*), y(*), z(*)
      real*8   xc(*), yc(*), zc(*)
      REAL*8   vx(*),vy(*),vz(*)
      REAL*8   FACEFLUX(3,*)
      REAL*8   Faces_ReferenceVector(3,*),Faces_NeumannFlux(*)
      REAL*8   BKCFLOWINT_NODE(*),ARENOD(*)
      REAL*8   FaceCentroid(3,*),AREFACE(*)
C     PASSED AND CALCULATED
      REAL*8   vn(nface)
C
C    e = face passing for (x1,y1,z1), (x2,y2,z2), (x3,y3,z3)
C    components of the normal to the face:
C    nx = alpha, ny=beta nz=gamma
C    alpha= (y2-y1)(z3-z1) - (z2-z1)(y3-y1)
C    beta = (x3-x1)(z2-z1) - (x2-x1)(z3-z1)
C    gamma= (x2-x1)(y3-y1) - (y2-y1)(x3-x1)
C    d= sqrt(alpha**2 + beta**2 + gamma**2)
C    normalized : nx= alpha/d, ny=beta/d, nz=gamma/d
C   
C   vn = integral along the face e of the inner product of v and n
C   vn = intg <v,n> ds = intg vx nx + vy ny + vz nz ds
C   setting z=f(x,y) the intg along the face e of a function 
C   g(x,y,z) can be computed as the integral into 2 variables x, y
C   intg g(x,y,f(x,y) sqrt(1 + dfx**2 + dfy**2) dx dy
C   In our case sqrt(1 + dfx**2 + dfy**2) = 1/|nz| and we have:
C   intg (vx nx + vy ny + vz nz ) 1/|nz| dx dy =
C   (vx nx + vy ny + vz nz ) 1/|nz| area(dx dy) = 
C   (vx nx + vy ny + vz nz ) 1/|nz| gamma/2 =
C   (vx nx + vy ny + vz nz ) d/2   =
C   (vx alpha + vy beta + vz gamma)/2
C        |          |          |
C       =          =           =
C   ( vx bi      vy ci       vz di )/2
C
C
       DO i=1,nface
          iel1= plist(1,i)
          IF(VELREC.NE.1)THEN
            do j=1,3
              xl(j) = x(iside(j,i))
              yl(j) = y(iside(j,i))
              zl(j) = z(iside(j,i))
            end do
            alpha = (yl(2) - yl(1))*(zl(3) - zl(1)) - 
     1              (yl(3) - yl(1))*(zl(2) - zl(1)) 
            beta  = (xl(3) - xl(1))*(zl(2) - zl(1)) - 
     1              (xl(2) - xl(1))*(zl(3) - zl(1)) 
            gamma = (xl(2) - xl(1))*(yl(3) - yl(1)) - 
     1              (xl(3) - xl(1))*(yl(2) - yl(1)) 
            prodscal =(xl(3)-xc(iel1))*alpha +(yl(3)-yc(iel1))*beta
     2              + (zl(3) - zc(iel1))*gamma
            if (prodscal.lt.0.d0) then
               alpha = -alpha
               beta = -beta
               gamma = -gamma
            end if
            iel2 = plist(2,i)
                 if (iel2.ne.0) then
                    vxface= 0.5d0*(vx(iel1) + vx(iel2))
                    vyface= 0.5d0*(vy(iel1) + vy(iel2))
                    vzface= 0.5d0*(vz(iel1) + vz(iel2))
                 else
                    vxface= vx(iel1)
                    vyface= vy(iel1)
                    vzface= vz(iel1)
                 end if
            vn(i) = 0.5d0*(vxface*alpha +vyface*beta +vzface*gamma)
          ELSE   
            FLUX(1)=FACEFLUX(1,i)
            FLUX(2)=FACEFLUX(2,i)
            FLUX(3)=FACEFLUX(3,i)
               DO J=1,4
                  IF(SIDE_CNC(J,iel1).EQ.i)THEN   
                     CALL EvalElementFaceNormal(iel1,J,
     1                    Normal1,Faces_ReferenceVector,
     2                    XC,YC,ZC,SIDE_CNC,FaceCentroid)
                  END IF
               END DO 
            vn(i)=DO3DSCALARPRODUCT(Normal1,FLUX)*AREFACE(I)
          END IF
        END DO
CC
CC CARLOTTA The following lines perform a correction of the boundary 
CC fluxes when velocity reconstruction is not performed. The 
CC fluxes are recovered starting from the back-calculated fluxes, 
CC after FLOW3D.   
CC         
        DO I=1,NFACE
           Node1=ISIDE(1,I)
           Node2=ISIDE(2,I)
           Node3=ISIDE(3,I)
           AREA(Node1)=0.0d0
           AREA(Node2)=0.0d0
           AREA(Node3)=0.0d0
           IF(PLIST(2,I).EQ.0)THEN
              IF(VELREC.NE.1)THEN
                  IF(Faces_FaceType(I).EQ.2)THEN 
                    VN(I)=Faces_NeumannFlux(I)*AREFACE(I)
                  END IF
                  IF(Faces_FaceType(I).EQ.0)THEN 
                    DO K=1,Node_FaceCount(Node1)
                        Faccia=Node_FaceIDS(Node1,K)
                          IF(PLIST(2,Faccia).EQ.0)THEN 
                            AREA(Node1)=AREA(Node1)+
     1                      AREFACE(Faccia)/3.0d0
                          END IF
                    END DO 
                    DO K=1,Node_FaceCount(Node2)
                        Faccia=Node_FaceIDS(Node2,K)
                          IF(PLIST(2,Faccia).EQ.0)THEN 
                            AREA(Node2)=AREA(Node2)+
     1                      AREFACE(Faccia)/3.0d0
                          END IF
                    END DO 
                    DO K=1,Node_FaceCount(Node3)
                        Faccia=Node_FaceIDS(Node3,K)
                          IF(PLIST(2,Faccia).EQ.0)THEN 
                            AREA(Node3)=AREA(Node3)+
     1                      AREFACE(Faccia)/3.0d0
                          END IF
                    END DO
                    VN(I)=-(BKCFLOWINT_NODE(Node1)/AREA(Node1)+
     1                     BKCFLOWINT_NODE(Node2)/AREA(Node2)+
     2                     BKCFLOWINT_NODE(Node3)/AREA(Node3))
     3                     *AREFACE(I)/3.0d0
                  END IF
              END IF   
           END IF 
        END DO
C
      RETURN
      END
