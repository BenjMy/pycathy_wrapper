      subroutine mlg_limiter(nnode,ntetra,nface,npfa_tra,n1,
     1                     flag_temp,neigh,
     2                     side_cnc,puntdir,plist,iside,ip4,
     3                     tetravert,tetra,
     3                     x,y,z,xc,yc,zc,conc_el,cdiff,presc,
     4                     precl,precr)
      implicit none
  
      INCLUDE 'CATHY.H' 
      integer ntetra,nnode, nface, npfa_tra,n1 
      integer flag_temp
      integer neigh(4,*),side_cnc(4,*)
      integer iside(3,*),plist(2,*)
      integer puntdir(*),ip4(4,4),tetravert(NMAX,*)
      integer tetra(4,*)

      integer i, j,l,itetra,itvd,ind(4),t,j1,j2,j3,p_lim
      integer ineigh,iface,idir,colnode,node,neartri
      
      integer idmax 
     
      real*8  x(*),y(*),z(*)
      real*8  xc(*),yc(*),zc(*)
      real*8  conc_el(*), cdiff(*), presc(*)
      real*8  ul(4), xl(4), yl(4), zl(4), xb,yb, zb, ubar,u2
      real*8  xm(4),ym(4),zm(4),xnode(3),ynode(3),znode(3),uxm(4)
      real*8  ux(5), uy(5),  uz(5), const(5) 
      real*8  x11,x21,x31,y12,y22,y32,z13,z23,z33
      real*8  u11,u12,u13, det_face, det, unorm(5),rm(4,5)
      real*8  ubtetra, lbtetra, poldn,lmax,lmin, alpha(4),alphalim(5)
      real*8  precl(*), precr(*), diff_u
      real*8  eps,zero
      real*8 der_x
      parameter (zero=0.d0,eps=1.e-10)

 
      logical tvd, noshoot, tbound, bordo

      do itetra = 1,ntetra
           bordo=.false.
C
C      xnode, ynode, znode: x, y and z-coordinates of nodes of 
C      each face
C
         do j = 1,4
            iface= side_cnc(j,itetra)
            do i=1,3
               xnode(i) = x(iside(i,iface))
               ynode(i) = y(iside(i,iface))
               znode(i) = z(iside(i,iface))
            end do
C
C   xm, ym, zm: coordinates of the centroid of each face
C
            call baric_face(xnode,ynode,znode,xm(j),ym(j),zm(j))
         end do

         do j = 1,4
C
C    assignment of values of conc_el, and x and z- coordinates of centroids
C    on the neighbouring triagles
C
            ineigh = neigh(j,itetra)
            if (ineigh.ne.0) then
                ul(j) = conc_el(ineigh) 
                xl(j) = xc(ineigh)
                yl(j) = yc(ineigh)
                zl(j) = zc(ineigh)
            else
C
C    on a boundary triangle
C
                iface = side_cnc(j,itetra)
                idir = puntdir(iface)
                if (idir.eq.0) then
                    ul(j) = conc_el(itetra)
C  CARLOTTA if it is not a Dirichlet boundary face, the 
C           concentration is 0
c                    ul(j) = 0.0d0 
                else
                    ul(j) = presc(idir)
                end if
                xl(j) = xm(j) 
                yl(j) = ym(j) 
                zl(j) = zm(j) 
            end if
         end do
C
C   values of conc_el and x and z-coordinates on itetra tetrehedron
C
         xb = xc(itetra)
         yb = yc(itetra)
         zb = zc(itetra)
         ubar = conc_el(itetra)
C
C   computation of gradients of 5  linear interpolation reconstruction
C

         do j = 1,4
            j1=ip4(j,1)
            j2=ip4(j,2)
            j3=ip4(j,3)
            x11= xl(j1) -xl(j2)
            x21= xl(j2) -xl(j3)
            x31= xl(j3) -xb
            y12= yl(j1) -yl(j2)
            y22= yl(j2) -yl(j3)
            y32= yl(j3) -yb
            z13= zl(j1) -zl(j2)
            z23= zl(j2) -zl(j3)
            z33= zl(j3) -zb
            u11= ul(j1) -ul(j2)
            u12= ul(j2) -ul(j3)
            u13= ul(j3) -ubar
            det_face= det(x11,y12,z13,x21,y22,z23,x31,y32,z33)
            ux(j)= det(u11,y12,z13,u12,y22,z23,u13,y32,z33)/det_face
            uy(j)= det(x11,u11,z13,x21,u12,z23,x31,u13,z33)/det_face
            uz(j)= det(x11,y12,u11,x21,y22,u12,x31,y32,u13)/det_face
            unorm(j) = dsqrt(ux(j)**2 +  uy(j)**2+ uz(j)**2)
            const(j) =ubar-ux(j)*xb-uy(j)*yb - uz(j)*zb
          end do
C  construction of the fifth interpolant, setting flag_interp=2 in the
C  iperplane subroutine
            call iperplane(2,xb,yb,zb,xl,yl,zl,ubar,ul,ux(5),uy(5),
     1                     uz(5),const(5))
            unorm(5) = dsqrt(ux(5)**2 +  uy(5)**2+ uz(5)**2)


            do l=1,5
             do j=1,4
             rm(j,l)=ux(l)*(xm(j)-xb)+uy(l)*(ym(j)-yb)+uz(l)*(zm(j)-zb)
               ineigh = neigh(j,itetra)
               if (ineigh.ne.0) then
                   u2 = conc_el(ineigh)
               else
                   iface = side_cnc(j,itetra)
                   idir = puntdir(iface)
                   if (idir.eq.0) then
C CARLOTTA if it is not a Dirichlet boundary face, the cocnentration is 0 
                  u2 = conc_el(itetra)
c                       u2 = 0.0d0 
                   else
                       u2 = presc(idir)
                   end if
               end if     
               diff_u= u2- ubar
            if (dabs(rm(j,l)).le.eps) rm(j,l)=zero
            if (dabs(u2-ubar).le.eps)   diff_u =zero
               if (rm(j,l).gt.max(diff_u,zero)) then
                  alpha(j) = max(diff_u,zero)/rm(j,l)
               else if (rm(j,l).lt.min(diff_u,zero)) then
                  alpha(j) = min(diff_u,zero)/rm(j,l)
               else
                  alpha(j) = 1.d0
               end if
             end do 
             alphalim(l) = min(alpha(1),alpha(2),alpha(3),alpha(4))   
             unorm(l)=alphalim(l)*dsqrt(ux(l)**2 + uy(l)**2+ uz(l)**2)  
            end do
            p_lim= idmax(5,unorm,1)
            do j=1,4
                uxm(j)= ubar + alphalim(p_lim)*rm(j,p_lim)
            end do  
             
               do j = 1,4
                 iface = side_cnc(j,itetra)
                 if(plist(1,iface).eq.itetra) then
                    precl(iface) = uxm(j)
                    if (flag_temp.eq.2) then
                      precl(iface)=precl(iface) +cdiff(itetra)
                    endif
c           
c on a boundary edge set also the `right' reconstructed value
c
                    if (plist(2,iface).eq.0) then
                        if(puntdir(iface).eq.0) then
C CARLOTTA If it is not a Dirichlet boundary face, the concentration is
C 0
C                             precr(iface) = precl(iface)
                               precr(iface) = 0.0d0
                             if (flag_temp.eq.2) then
c                               precr(iface)=precr(iface)+cdiff(itetra)
                               precr(iface)=0.0d0
                             endif
                        else
                             precr(iface) = presc(puntdir(iface))
                        end if
                    end if
                 else
                    precr(iface) = uxm(j)
                    if (flag_temp.eq.2) then
                      precr(iface)=precr(iface) +cdiff(itetra)
                    endif
                 end if
              end do           
      end do
 
 
 
      return
      end                      
