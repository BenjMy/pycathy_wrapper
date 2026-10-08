C  *********************** IPERPLANE *************************************
C  construction of the surface L(x,y,z) necessary to obtain second order
C  accuracy  in the FV reconstruction
C ************************************************************************
C
      subroutine iperplane(flag_shoot, xb,yb,zb,xl,yl,zl,
     1                     ubar,ul, ux,uy,uz,const)

      implicit none
      integer flag_shoot, j, info, ipiv(4)
      real*8  xb,yb,zb, xl(4), yl(4), zl(4), ubar, ul(4)
      real*8  ux,uy,uz,const
      real*8  x11,x21,x31,y12,y22,y32,z13,z23,z33,a(4,4),b(4)
      real*8  u11,u12,u13, det_face, det


      if (flag_shoot.eq.2) then
C  interpolant constructed considering the 4 neighbouring tetrahedra of the
C  refering tetrahedron
            x11= xl(1) -xl(2)
            x21= xl(2) -xl(3)
            x31= xl(3) -xl(4)
            y12= yl(1) -yl(2)
            y22= yl(2) -yl(3)
            y32= yl(3) -yl(4)
            z13= zl(1) -zl(2)
            z23= zl(2) -zl(3)
            z33= zl(3) -zl(4)
            u11= ul(1) -ul(2)
            u12= ul(2) -ul(3)
            u13= ul(3) -ul(4)
            det_face= det(x11,y12,z13,x21,y22,z23,x31,y32,z33)
            ux= det(u11,y12,z13,u12,y22,z23,u13,y32,z33)/det_face
            uy= det(x11,u11,z13,x21,u12,z23,x31,u13,z33)/det_face
            uz= det(x11,y12,u11,x21,y22,u12,x31,y32,u13)/det_face
            const =ubar-ux*xb-uy*yb - uz*zb
       else if (flag_shoot.eq.4) then
C  MLS interpolant: 
C  interpolant minimizing the values at the 4 neighbours and at the refering
C  tetrahedron with the least square method
            a(1,1)= xb**2
            a(2,1)= xb*yb
            a(3,1)= xb*zb
            a(2,2)= yb**2
            a(3,2)= yb*zb
            a(3,3)= zb**2
            a(1,4)= xb
            a(2,4)= yb
            a(3,4)= zb
            a(4,4) = 5.0d0
            b(1) = ubar*xb
            b(2) = ubar*yb
            b(3)= ubar*zb
            b(4)=ubar
            do j=1,4
             a(1,1)=a(1,1) + xl(j)**2 
             a(2,1) = a(2,1) + xl(j)*yl(j)
             a(3,1)= a(3,1) + xl(j)*zl(j)
             a(2,2)= a(2,2) + yl(j)**2
             a(3,2)= a(3,2) + yl(j)*zl(j)
             a(3,3) =  a(3,3)+ zl(j)**2
             a(1,4)= a(1,4)+xl(j)
             a(2,4)= a(2,4)+yl(j)
             a(3,4)= a(3,4)+zl(j)
             b(1)= b(1) + ul(j)*xl(j)
             b(2)= b(2) + ul(j)*yl(j)
             b(3)= b(3) + ul(j)*zl(j)
             b(4) = b(4) + ul(j)
            end do
             a(1,2) = a(2,1)
             a(1,3)=a(3,1)
             a(2,3)=a(3,2)
             a(4,1)=a(1,4)
             a(4,2)=a(2,4)
             a(4,3)=a(3,4)
            call  dgetrf(4,4,a,4,ipiv,info)
            call  dgetrs('N',4,1,a,4,ipiv,b,4,info) 
            ux=b(1)      
            uy=b(2)
            uz=b(3)
            const =ubar-ux*xb-uy*yb - uz*zb
           end if

         return
         end                 
