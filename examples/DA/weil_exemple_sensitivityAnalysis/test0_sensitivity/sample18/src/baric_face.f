C
C**************************  BARIC_FACE **********************************
C
C  calculate the coordinates of the centroid of a face of a tetrahedron
C
C***********************************************************************
C
      subroutine baric_face(xl,yl,zl,xb,yb,zb)
 
      IMPLICIT  NONE

      REAL*8    xl(3), yl(3), zl(3),xb, yb, zb
 
      xb = (xl(1) + xl(2) + xl(3))/3.d0
      yb = (yl(1) + yl(2) + yl(3))/3.d0
      zb = (zl(1) + zl(2) + zl(3))/3.d0
 
      RETURN
      END
