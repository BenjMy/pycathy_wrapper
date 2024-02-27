C
C*************************** nodetotetra *******************************
C
C  FOR TRANSPORT: pass the values from the centroid of each triangle to
C  the nodes of the triangle
C
C***********************************************************************
C
      SUBROUTINE nodetotetra(n,nt,tetra,pold,pnold,volu,volnod,
     1                       pnodi,pel)
C
      IMPLICIT  NONE
      INTEGER   itria,j,jj,ii,k,i,inod
      INTEGER   n, nt, tetra(5,*)
      REAL*8    pold(*), pnold(*)
      REAL*8    volu(*),volnod(*),pel(*),pnodi(*)
      real*8    volneigh(n),cont(n)
C
      call init0r(nt,pold)
      call init0r(n,cont)
      call init0r(n,volneigh)
c      do itria = 1, nt
c          do j=1,4
c             jj=tetra(j,itria)
c             pold(itria) = pold(itria) + pnold(jj)
c          end do
c      end do
c      do j=1,nt
c         pold(j) = pold(j)/4.0d0
c      end do
         
      DO K=1,NT
         DO II=1,4
            INOD=TETRA(II,K)
            volneigh(INOD)=volneigh(INOD)+VOLU(k)*pel(k)
            cont(inod)=cont(inod)+1
         END DO
      END DO
         
        do i = 1, nt
        DO II=1,4 
        inod=tetra(ii,i)
        pold(i) = pold(i) + pnold(inod)*volnod(inod)*pnodi(inod)
     1                    *volu(i)*pel(i)/volneigh(inod)
        end do
         pold(i) = pold(i)/ (volu(i) * pel(i))
      end do
C
      RETURN
      END
