C  ***************************  HERON  *************************
C  function that computes the area of faces of each element 
C   by applying  Heron's formula:
C   given a,b,c the length of the edge of a triangle, 
C   set s=(a+b+c)/2 the formula is:
C    area = dsqrt(s*(s-a)*(s-b)*(s-c))
C
C **************************************************************
       real*8 function heron(xfac,yfac,zfac)
       implicit none
       integer j
       real*8 xfac(3), yfac(3), zfac(3), length(3)
       real*8 s

       do j=1,3
         length(j) =dsqrt( (xfac(j)- xfac(mod(j,3)+1) )**2 +
     1              (yfac(j) - yfac(mod(j,3)+1) )**2 +
     2              (zfac(j) - zfac(mod(j,3)+1) )**2  )
       end do
       s= ( length(1) + length(2) +length(3) )*0.5d0
       heron = s
       do j=1,3
         heron =heron*(s-length(j))
       end do
       heron = dsqrt(heron)
       return
       end
         

