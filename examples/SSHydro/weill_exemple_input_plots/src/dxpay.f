      subroutine dxpay(n,dx,incx,da,dy,incy)
c
c     constant times a vector plus a vector.
c     y:=x+ay
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      implicit none
      integer i,incx,incy,ix,iy,m,mp1,n
cxcx  real*8 dx(n),dy(n),da
      real*8 dx(*),dy(*),da
c
      if(n.le.0)return
      if (da .eq. 0.0d0) then
         call dcopy(n,dx,1,dy,1)
         return
      end if
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(iy) + da*dy(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i) + da*dy(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dx(i) + da*dy(i)
        dy(i + 1) = dx(i + 1) + da*dy(i + 1)
        dy(i + 2) = dx(i + 2) + da*dy(i + 2)
        dy(i + 3) = dx(i + 3) + da*dy(i + 3)
   50 continue
      return
      end
