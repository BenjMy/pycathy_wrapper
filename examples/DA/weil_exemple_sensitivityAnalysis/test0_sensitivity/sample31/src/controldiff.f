C    ********************* CONTROLDIFF *******************************  
C    to control if there is diffusion on the transport equation
C
C    *****************************************************************
C
C
      subroutine controldiff(nstr,nzone,alfal,alfat,diffus,cdiffus)
C
C     cdiffus  is a parameter to check if there is diffusion or not
C     cdiffus > 1 ---> diffusion
C     cdiffus = 0 ---> no diffusion
C
      implicit none
      integer nstr,nzone, cdiffus
      integer i,j
      real*8  alfal(nstr,nzone), alfat(nstr,nzone), diffus

      cdiffus = 0
      do i=1,nstr
          do j=1,nzone
            if (alfal(i,j).ne.0.0)  cdiffus = cdiffus +1
            if (alfat(i,j).ne.0.0)  cdiffus = cdiffus +1
          end do
      end do

      if (diffus.ne.0.0) cdiffus = cdiffus + 1

      return
      end
