C
C************************ NORMRES ***********************************
C
      real*8 function normres(nequ,ndir,noddir,vec,scr)
c
c  calcultes the Euclidian length of vector vec, 
c  not considering the "ndir" components given in vector "noddir"
c
      implicit none
      integer  nequ,ndir
      integer  i
cxcx  integer  noddir(ndir)
      integer  noddir(*)
cxcx  real*8   vec(nequ),scr(nequ)
      real*8   vec(*),scr(*)
      real*8   zero
      real*8   dnrm2

      parameter (zero=0.d0)

      if (ndir.eq.0) then
         normres=dnrm2(nequ,vec,1)
      else
         call dcopy(nequ,vec,1,scr,1)
         do i=1,ndir
            scr(noddir(i))=zero
         end do
         normres=dnrm2(nequ,scr,1)
      end if

      return
      end
