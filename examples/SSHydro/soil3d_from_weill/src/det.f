c
c    ************************  DET *******************************
c
c              calculate the determinant of matrix A(3,3)
c     *************************************************************
       real*8 function det(a11,a12,a13,a21,a22,a23,a31,a32,a33)
       implicit none
       real*8 a11,a12,a13,a21,a22,a23,a31,a32,a33
       det=a11*(a22*a33-a23*a32)+a12*(a23*a31-a21*a33)+
     1     a13*(a21*a32-a22*a31)
       return
      end
