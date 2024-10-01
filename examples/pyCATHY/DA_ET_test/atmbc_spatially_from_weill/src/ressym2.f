C
C**************************  RESSYM ************************************
C
C  calculate the residual from the conjugate gradients
C  solution of a symmetric linear system
C
C***********************************************************************
C
      subroutine ressym2(nequ,nterm,ndir,ia,ja,noddir,sysmat,rhs,xk,res)

      implicit  none
      integer   i,k,m,mm
      integer   nequ,nterm,ndir
cxcx  integer   ia(nequ+1),ja(nterm),noddir(ndir)
      integer   ia(*),ja(*),noddir(*)

cxcx  real*8    sysmat(nterm),rhs(nequ),xk(nequ),res(nequ)
      real*8    sysmat(*),rhs(*),xk(*),res(*)
c
      do k=1,nequ
         res(k)=rhs(k)
      end do
cxcx  do k=1,ndir
cxcx     res(noddir(k))=0.0d0
cxcx  end do
      if (ndir .gt. 0) then
         do k=1,ndir
            res(noddir(k))=0.0d0
         end do
      end if
      do k=1,nequ
         m=ia(k)
         mm=ia(k+1)-1
         res(k)=res(k)-sysmat(m)*xk(k)
         do i=m+1,mm
            res(k)=res(k)-sysmat(i)*xk(ja(i))
            res(ja(i))=res(ja(i))-sysmat(i)*xk(k)
         end do
      end do
c
      return
      end
