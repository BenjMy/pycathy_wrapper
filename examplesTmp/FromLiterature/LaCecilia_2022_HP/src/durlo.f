C
C**************************  DURLO  ***********************************
C
C     3D  advective flux for the FINITE VOLUME scheme
C***********************************************************************
C
      subroutine durlo(ntetra,nface,volur,plist,vn,
     1                precl,precr,time,flux)

      implicit none

      INCLUDE 'CATHY.H'
      integer itetra,iface,itetra1,itetra2
      integer ntetra,nface
      integer plist(2,*)
      real*8 cup,qn
      real*8 time
cxcx  real*8 volur(*),c,a
      real*8 volur(*)
      real*8 vn(*),precl(*),precr(*),flux(*)
C
      call init0r(ntemax,flux)
c
cxcx "c" computed below is never used ...
cxcx   do iface=1,nface
cxcx      c=c+precl(iface)
cxcx   end do
cxcx   do iface=1,nface
cxcx      c=c+precr(iface)
cxcx   end do
      do iface = 1,nface
         itetra1 = plist(1,iface)
         itetra2 = plist(2,iface)
          qn=vn(iface)
c     solve the Riemann problem
c
         if (qn.gt.0.0) then
            cup = precl(iface)
         else
            cup = precr(iface)
         end if
         flux(itetra1)=flux(itetra1)+qn*cup
         if(itetra2.ne.0) flux(itetra2)=flux(itetra2)-qn*cup
      end do
      do itetra = 1,ntetra
         flux(itetra)= -dabs(volur(itetra))*flux(itetra)
      end do
C
      return
      end
