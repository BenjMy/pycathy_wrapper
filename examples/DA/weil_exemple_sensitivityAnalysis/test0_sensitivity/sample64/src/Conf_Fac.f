      logical function conf_fac(face1,face2)
C-----------------------------------------------------------------------
C  Function che restituisce il valore vero se face1 < face2
C-----------------------------------------------------------------------
      integer  face1(3),face2(3)
C
      conf_fac = .True.
      if (face1(1) .gt. face2(1)) then
         conf_fac = .false.
      elseif (face1(1) .eq. face2(1)) then
         if (face1(2) .gt. face2(2)) then
            conf_fac = .false.
         elseif (face1(2) .eq. face2(2)) then
            conf_fac = (face1(3) .le. face2(3))
         endif
      endif
C
      return
C
      end

