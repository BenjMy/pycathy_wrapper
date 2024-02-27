      real*8 function norma2(n,x)
      implicit none
      integer n,i
      real*8  x(n)

      norma2=0.d0
      do i=1,n 
          norma2=norma2+ x(i)*x(i)
      end do
      norma2=dsqrt(norma2)

      return
      end
