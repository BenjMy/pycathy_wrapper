      subroutine Swap_t(a,b)
c
      implicit none
c
      integer  a(4),b(4)
      integer  i,temp
c
      do i=1,4
         temp = a(i)
         a(i) = b(i)
         b(i) = temp
      end do 
c
      return
c
      end
