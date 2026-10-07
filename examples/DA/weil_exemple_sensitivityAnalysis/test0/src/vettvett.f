      real*8 function vettvett(nequ,x,y)
      implicit none
      integer nequ,i
      real*8 x(nequ), y(nequ)
      vettvett=0.d0
      do i=1,nequ
        vettvett=vettvett+x(i)*y(i)
      end do
 
      return
      end
