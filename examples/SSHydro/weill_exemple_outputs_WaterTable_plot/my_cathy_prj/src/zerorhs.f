C
C************************* ZERORHS *************************************
C
      subroutine zerorhs(nequ,info,ndir,bnorm,resnorm,sol,rhs,exit_test)
 
      implicit none

      integer   nequ,info,ndir,i

      real*8    bnorm,resnorm,zero
cxcx  real*8    sol(nequ),rhs(nequ)
      real*8    sol(*),rhs(*)
      real*8    dnrm2

      parameter (zero = 0.d0)

      logical   exit_test

c
c handles case of zero RHS
c
      exit_test=.FALSE.
      if (bnorm.eq.zero) then
         if (ndir.eq.0) then
            do i = 1,nequ
               sol(i) = zero
            end do
            resnorm = zero
            info = 0
            exit_test=.true.
         else if (dnrm2(nequ,rhs,1).gt.zero) then
            bnorm = float(nequ)
         else
            do i = 1,nequ
               sol(i) = zero
            end do
            resnorm = zero
            info = 0
            exit_test=.true.
         end if
      end if

      return 
      end
