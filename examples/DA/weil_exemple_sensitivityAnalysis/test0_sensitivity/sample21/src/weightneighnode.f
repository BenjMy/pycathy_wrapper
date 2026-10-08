C
C**********************  weightneighnode *******************************
C
C  calculate for each node the value
C  1/r_i where r_i is the distance of the node to the centroid
C  of element T_i, where T_i is a neighboring element of the node
C  this distancs will be used as  a weighting
C
C***********************************************************************
C
      SUBROUTINE weightneighnode(N,NT,x,y,z,xbar,ybar,zbar,TETRA,
     1                          disnod,weightnod)
C
      IMPLICIT  NONE
      INTEGER   i,j,k,elle,inod
      INTEGER   N,NT
      INTEGER   TETRA(5,*),cont(N)
      real*8    x(*),y(*),z(*),xbar(*),ybar(*),zbar(*)
      real*8    disnod(24,*),weightnod(*),aux
C
      do i=1,n
        cont(i)=1
        weightnod(i)=0.d0
        do j=1,24
         disnod(j,i)=0.d0
        end do
      end do
c     write(118,*) 'nodo', 'tetraedro'
      DO k=1,NT
         DO elle=1,4
           inod=TETRA(elle,K)
           aux=dsqrt((x(inod)-xbar(k))**2+(y(inod)-ybar(k))**2 +
     1            (z(inod)-zbar(k))**2 )
           aux=1.0d0/aux
           disnod(cont(inod),inod)=aux
           cont(inod)=cont(inod)+1
           weightnod(inod)=weightnod(inod)+aux
         END DO
      END DO


      RETURN
      END
