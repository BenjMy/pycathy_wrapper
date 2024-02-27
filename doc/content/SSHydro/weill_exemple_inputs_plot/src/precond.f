C
C************************ PRECOND ***********************************
C
      subroutine precond(iprec,nequ,ntermp,iap,jap,prec,
     1                   vec,scr,pvec)
c
c applies the preconditioner to a vector
c
c    pvec:=prec*vec
c
c  iprec = 0 -> no preconditioning
c          1 -> diagonal preconditioning
c          2 -> modified diagonal preconditioning
c          3 -> Choleski preconditioning
c          4 -> AINV preconditioning

      implicit none
      integer  iprec,nequ,ntermp
      integer  i
cxcx  integer  iap(nequ+1),jap(ntermp)
      integer  iap(*),jap(*)

cxcx  real*8   prec(ntermp),vec(nequ),pvec(nequ)
cxcx  real*8   scr(nequ)
      real*8   prec(*),vec(*),pvec(*)
      real*8   scr(*)

      if (iprec.eq.0) then
         call dcopy(nequ,vec,1,pvec,1)
         return
      else if(iprec.eq.1 .or. iprec.eq.2) then
         do i=1,nequ
            pvec(i)=prec(i)*vec(i)
         end do
         return
      else if(iprec.eq.3) then
         call lsolve(nequ,ntermp,iap,jap,prec,vec,pvec)
         return
      else if(iprec.eq.4) then
         call ztzvec(nequ,ntermp,iap,jap,prec,vec,scr,pvec)
         return
      else 
         write(*,*) ' wrong iprec value',iprec
         stop
      end if

      end
