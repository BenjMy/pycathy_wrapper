      subroutine tvd_durlo(nnode,ntetra,nface,npfa_tra,n1,
     1                     flag_limiter,flag_temp,neigh,
     2                     side_cnc,puntdir,plist,iside,ip4,
     3                     tetravert,tetra,
     3                     x,y,z,xc,yc,zc,conc_el,cdiff,
     4                     prescfa_tra,precl,precr)
      implicit none
  
      INCLUDE 'CATHY.H' 
      integer flag_limiter,ntetra,nnode,nface,npfa_tra,n1 
      integer flag_temp
      integer neigh(4,*),side_cnc(4,*)
      integer iside(3,*),plist(2,*)
      integer puntdir(*),ip4(4,4),tetravert(nmax,*)
      integer tetra(4,*)
      integer i,j,itetra,itvd,ind(4),t,j1,j2,j3
      integer ineigh,iface,idir,colnode,node,neartri
      real*8  x(*),y(*),z(*)
      real*8  xc(*),yc(*),zc(*)
      real*8  conc_el(*),cdiff(*),prescfa_tra(*)
      real*8  ul(4),xl(4),yl(4),zl(4),xb,yb,zb,ubar,u2
      real*8  xm(4),ym(4),zm(4),xnode(3),ynode(3),znode(3),uxm(4)
      real*8  ux(4),uy(4),uz(4),const(4) 
      real*8  x11,x21,x31,y12,y22,y32,z13,z23,z33
      real*8  u11,u12,u13,det_face,det,unorm(4)
      real*8  ubtetra,lbtetra,poldn,lmax,lmin
      real*8  precl(*),precr(*)  
      real*8  eps
      real*8  der_x
      parameter (eps=1.e-10)
 
      logical tvd, noshoot, tbound
      
      do itetra = 1,ntetra
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
c  If it is not a boundary Dirichlet face, then the 
c  concentration is 0
c                     ul(j) = 0.0d0 
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
C   computation of gradients of the three linear interpolation reconstruction
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
            ind(j) = j
          end do
C
C    computation of the maximum gradient
C    ind(1)---> index of maximum gradient
C    ind(2)---> index of second largest gradient
C    ind(3)---> index of third  largest gradient
C    ind(4)---> index of minimum gradient
C
       if (unorm(ind(1)).lt.unorm(ind(2))) call swap(ind,1,2)
       if (unorm(ind(1)).lt.unorm(ind(3))) call swap(ind,1,3)
       if (unorm(ind(1)).lt.unorm(ind(4))) call swap(ind,1,4)
       if (unorm(ind(2)).lt.unorm(ind(3))) call swap(ind,2,3)
       if (unorm(ind(2)).lt.unorm(ind(4))) call swap(ind,2,4)
       if (unorm(ind(3)).lt.unorm(ind(4))) call swap(ind,3,4)
C
C      we check if the linear interpolant with maximum gradient calculated 
C      on the centroid of each face is between the values of conc_el on the 
C      itetra tetrahedron and the value of conc_el on the adiacent tetrahedra.
C      If this is not true for the linear interpolant with maximum gradient
C      we check a control for the second largest to the minimum
C      
      if (flag_limiter.eq.3) then
C averaged gradient
C no limiter will be  considered
           ux(ind(1))=(ux(ind(1))+ux(ind(2))+
     1                 ux(ind(3))+ux(ind(4)))/4.0d0
           uy(ind(1))=(uy(ind(1))+uy(ind(2))+
     1                 uy(ind(3))+uy(ind(4)))/4.0d0
           uz(ind(1))=(uz(ind(1))+uz(ind(2))+
     1                 uz(ind(3))+uz(ind(4)))/4.0d0
           const(ind(1)) =ubar-ux(ind(1))*xb -uy(ind(1))*yb-
     1                   uz(ind(1))*zb
          do j=1,4 
                uxm(j) = ux(ind(1))*xm(j) + uy(ind(1))*ym(j)+
     1                   uz(ind(1))*zm(j) + const(ind(1))
          end do 
C
      elseif (flag_limiter.eq.2) then
C maximum gradient
C no limiter will be considered
          do j=1,4 
                uxm(j) = ux(ind(1))*xm(j) + uy(ind(1))*ym(j)+
     1                   uz(ind(1))*zm(j) + const(ind(1))
          end do 
C
      elseif (flag_limiter.eq.1 ) then 
C
C     see Durlofsky, Engquist, Osher, JCP 98, 64-73 (1992)
C     on page  66, second column, after equation (2.7)
C     ... By analogy with limiting procedure on one space dimension
C     a valid non-compressive limiter corresponds to the selection of
C     the L_i for which unorm is the minimum. ...
C     At the extrema, a first-order approximation is used.
C     This limiter is applicable to problems involving nonlinear flux function.
C    
C    we use the limiter for one-dimensional problems
          do j=1,4 
                uxm(j) = ux(ind(4))*xm(j) + uy(ind(4))*ym(j)+
     1                   uz(ind(4))*zm(j) + const(ind(4))
           end do
           if (dabs(ubar).gt.max(dabs(ul(1)),dabs(ul(2)),
     1                           dabs(ul(3)),dabs(ul(4)))) 
     1     then
              do j=1,4
                 uxm(j)=ubar
              end do
           end if
      else
C
C  Durlofsky's procedure extended to 3d
C
           tvd = .false.
           i = 1
           do while(.not.tvd.and.i.le.4)
             itvd = 0
             j=0
             do while (j.lt.4.and.itvd.eq.j)
                j=j+1
                uxm(j) = ux(ind(i))*xm(j) + uy(ind(i))*ym(j) +
     1                   uz(ind(i))*zm(j) + const(ind(i))
                ineigh = neigh(j,itetra)
                if (ineigh.ne.0) then
                    u2 = conc_el(ineigh)
                else
                    iface = side_cnc(j,itetra)
                    idir = puntdir(iface)
                    if (idir.eq.0) then
                        u2 = conc_el(itetra)
c if it is not a Dirichlet boundary face, the concentration is 0
c                         u2 = 0.0d0
                    else
                        u2 = prescfa_tra(idir)
                    end if
                end if
                 if ( (uxm(j)-u2)*(uxm(j)-ubar).le.eps) itvd=itvd+1
             end do
             tvd = (itvd.eq.4)
             i = i+1 
           end do
C
        end if
C
C     if tvd is true, i.e. the Durlofsky reconstruction is sufficient or
C     if noshoot is true, i.e. the Liu reconstruction goes well, then  we
C     use the linear interpolant to obtain the new value on each tetrahedron.
C
C     if the logical variable tvd is true, this means that one of the three 
C     linear reconstruction  satisfies the requirements of tvd
C
C     if tvd is true, we use the linear interpolant to obtain the
C     new value on each triangle.
C 
C     in precl and in precr we put the values of the linear reconstruction
C     on the midpoint of iface, inside and outside triang, 
C     by following  the normal direction indicated  by plist 
C
         if (flag_limiter.ne.0.or.tvd) then
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
c                             precr(iface) = precl(iface)
c  CARLOTTA if it is not a boundary Dirichlet face, the concentration is 0  
                              precr(iface) = 0.0d0 
                            if (flag_temp.eq.2) then
                              precr(iface) = 0.0d0 
c                               precr(iface)=precr(iface)+cdiff(itetra)
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
         else
C
C   if tvd is not true, then we operate by a constant
C   reconstruction as the first-order godunov procedure
C questo si ha con LIu ma con la procedura di Durlofsky
C
              do j = 1,4
                 iface = side_cnc(j,itetra)
                 if(plist(1,iface).eq.itetra) then
                    precl(iface) = ubar 
                    if (flag_temp.eq.2) then
                      precl(iface)=precl(iface) +cdiff(itetra)
                    end if
 
c on a boundary edge we set also the `right' reconstructed value
c
                    if (plist(2,iface).eq.0) then
                        if(puntdir(iface).eq.0) then
c                             precr(iface) = ubar
c CARLOTTA if it is not a boundary DIrichelt face, the right value is 0 
                             precr(iface) = 0.0d0
                            if (flag_temp.eq.2) then
c                            precr(iface)=precr(iface) +cdiff(itetra)
                            precr(iface)=0.0d0
                            end if
                        else
                             precr(iface) = prescfa_tra(puntdir(iface))
                        end if
                    end if
                 else
                    precr(iface) = ubar 
                    if (flag_temp.eq.2) then
                       precr(iface)=precr(iface) +cdiff(itetra)
                    end if
                 end if
              end do
          end if
      end do

      return
      end

C
C
C

      subroutine swap(vec,i,j)
     
      implicit none
      integer i,j,vec(4),temp

      temp = vec(i)
      vec(i) = vec(j) 
      vec(j) = temp  

      return
      end
