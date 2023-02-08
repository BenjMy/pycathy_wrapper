      subroutine extr_tvd(nnode,ntetra,nface,npfa_tra,n1,
     1                     flag_limiter,flag_interp,flag_temp,neigh,
     2                     side_cnc,puntdir,plist,iside,ip4,
     3                     tetravert,tetra,
     3                     x,y,z,xc,yc,zc,conc_el,cdiff,prescfa_tra,
     4                     precl,precr)
      implicit none
   
      INCLUDE 'CATHY.H'
      integer flag_limiter,ntetra,nnode,nface,npfa_tra,n1 
      integer flag_interp,flag_temp
      integer neigh(4,*),side_cnc(4,*)
      integer iside(3,*),plist(2,*)
      integer puntdir(*), ip4(4,4),tetravert(NMAX,*)
      integer tetra(4,*)

      integer i, j,itetra,itvd,ind(4),t,j1,j2,j3
      integer ineigh, iface,idir,colnode,node,neartri

      real*8  x(*),y(*),z(*)
      real*8  xc(*),yc(*),zc(*)
      real*8  conc_el(*),prescfa_tra(*)
      real*8  cdiff(*)
      real*8  ul(4),xl(4),yl(4),zl(4),xb,yb,zb,ubar,u2
      real*8  xm(4),ym(4),zm(4),xnode(3),ynode(3),znode(3),uxm(4)
      real*8  ux,uy,uz,const,rm(4)
      real*8  x11,x21,x31,y12,y22,y32,z13,z23,z33,a(4,4),b(4),xsol(4)
      real*8  u11,u12,u13,det_face,det,unorm,umax,umin
      real*8  ubtetra,lbtetra,poldn,lmax,lmin
      real*8  precl(*),precr(*),alpha(4),alphalim,diff_u
      real*8  eps,zero,c,e,f,g
      real*8  der_x
      parameter (eps=1.e-3,zero=0.d0)

 
      logical tvd, noshoot, tbound

     
       do itetra = 1,ntetra
C
C      xnode, ynode, znode: x, y and z-coordinates of nodes of 
C      each face

         do j = 1,4
            iface= side_cnc(j,itetra)
            do i=1,3
               xnode(i) = x(iside(i,iface))
               ynode(i) = y(iside(i,iface))
               znode(i) = z(iside(i,iface))
            end do

C   xm, ym, zm: coordinates of the centroid of each face
C

           call baric_face(xnode,ynode,znode,xm(j),ym(j),zm(j))
         end do

         do j = 1,4
C
C    assignment of values of conc_el, and x and z- coordinates of centroids
C    on the neighbouring element
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
C CARLOTTA if it is not a boundary Dirichlet face, then the right 
c          concentration is 0
                    ul(j) = conc_el(itetra)
c                    ul(j) = 0.0d0 
                else
                    ul(j) = prescfa_tra(idir)
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
C   computation of gradients of the linear interpolation reconstruction
C   taking into account the four neighbours of itetra (flag_interp.eq.2)
C   or considering the method of least square (flag_interp.eq.4)
C
         call iperplane(flag_interp,xb,yb,zb,xl,yl,zl,ubar,ul,
     1                  ux,uy,uz,const)
 
           
             do j=1,4
                rm(j)=ux*(xm(j)-xb)+uy*(ym(j)-yb)+uz*(zm(j)-zb)
                uxm(j) = ux*xm(j) + uy*ym(j)+ uz*zm(j) + const
                ineigh = neigh(j,itetra)
                if (ineigh.ne.0) then
                    u2 = conc_el(ineigh)
                else
                    iface = side_cnc(j,itetra)
                    idir = puntdir(iface)
                    if (idir.eq.0) then
                        u2 = conc_el(itetra)
C CARLOTTA if it is not a Dirichlet boundary face, the right 
C concentration is 0
c                         u2 = 0.0d0 
                    else
                         u2 = prescfa_tra(idir)
                    end if
                end if
               if (flag_limiter.eq.1) then
C  extremum limiter (simple form)
                 if (uxm(j).gt.max(ubar,u2)) then
                      uxm(j) = max(ubar,u2)
                 elseif (uxm(j).lt.min(ubar,u2) ) then
                      uxm(j) = min(ubar,u2)
                 end if
               else if (flag_limiter.eq.2) then
               if (dabs(rm(j)).le.eps) rm(j)=zero
               diff_u=u2-ubar
               if (dabs(diff_u).le.eps) diff_u=zero
                 if (rm(j).gt.max(diff_u,zero)) then
                    alpha(j) = max(diff_u,zero)/rm(j)
                 else if (rm(j).lt.min(diff_u,zero)) then
                    alpha(j) = min(diff_u,zero)/rm(j)
                 else
                    alpha(j)=1.0d0
                 end if
               else if (flag_limiter.eq.4) then

               if (dabs(rm(j)).le.eps) rm(j)=zero
               diff_u=uxm(j) - ubar
               umax= max(ubar, ul(1), ul(2), ul(3), ul(4))
               umin= min(ubar, ul(1), ul(2), ul(3), ul(4))
               if (dabs(diff_u).le.eps) diff_u=zero
                 if (diff_u.eq.zero) then
                  alpha(j)=1.0d0
                 else if (diff_u.gt.zero) then
                    alpha(j) = min( 1.0d0, (umax-ubar)/rm(j))
                 else if (diff_u.lt.zero) then
                    alpha(j) = min(1.0d0, (umin-ubar)/rm(j))
                  end if
c              write(1115,*)(umin-ubar),(umin-ubar)/rm(j),alpha(j)
               end if 
             end do
            if ((flag_limiter.eq.2).or.(flag_limiter.eq.4)) then
              alphalim=min(alpha(1), alpha(2), alpha(3), alpha(4))
              do j=1,4
                uxm(j) = ubar + alphalim*rm(j)
              end do
            else if (flag_limiter.eq.3) then
              alphalim=1.0d0
              do j=1,4
                uxm(j) = ubar + alphalim*rm(j)
              end do
            end if 
C
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
C  CARLOTTA If it is not a Dirichelt boundary face, the right
C  concentration is 0
c                             precr(iface) = precl(iface)
                             precr(iface) = 0.0d0 
                             if (flag_temp.eq.2) then
c                               precr(iface)=precr(iface)+cdiff(itetra)
                               precr(iface)=0.0d0
                             endif
                        else
                             precr(iface) = prescfa_tra(puntdir(iface))
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
       
c         write(1115,*)'UL',C

     

      return
      end
