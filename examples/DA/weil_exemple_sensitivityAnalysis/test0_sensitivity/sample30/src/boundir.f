C
C************************** BOUNDIR  ************************************
C
C puntdir are the pointers to pass from a global dirichlet numbering
C to a local dirichlet numbering
C***********************************************************************
C
      subroutine BOUNDIR(nface,np,contp,puntdir)
      implicit none
      integer nface, np
      integer contp(*), puntdir(*)
      integer i,iface
      
C     puntdir is the inverse function of edgedir
C
      do i=1,nface
         puntdir(i)=0
      end do
      do i=1,np
         iface=contp(i)
         puntdir(iface)=i
      end do

      return
      end
