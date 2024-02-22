c
c************************* ISTRATO *****************************
c
c  individua lo strato cui appartiene un determinato tetraedro
c
c*******************************************************
      integer function istrato(itetra,ntria)
      implicit none

      integer itetra,ntria,ir
      istrato=1+itetra/(ntria*3)
      ir=mod(itetra,ntria*3)
      if(ir.eq.0) then
       istrato=istrato-1
      endif
      return
      end
