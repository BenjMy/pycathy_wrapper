C
C************************* PCG **************************************
C
      subroutine pcg(iprt,iexit,iprec,imax,nequ,nterm,ntermp,ndir,
     1               info,iter,noddir,ia,ja,iap,jap,tolcg,bnorm,
     2               resini,resnorm,sysmat,prec,rhs,sol,scr,res,pk)
c
c  preconditioned conjugate gradient method
c  (cfr. Gambolati 1994, page 148)
c
c  uses the symmetric CSR storage scheme for the sparse matrix
c
c  Integer scalars
c  ***************
c
c iprt         = 0 does not print the convergence profile 
c                  (n. of iteration, PCG relative residual norm)
c              > 0 prints the convergence profile in unit "iprt"
c iprec        = 0   no preconditioning
c                1   D^-1 preconditioning
c                2   D^-1-norm preconditioning
c                3   Choleski preconditioning
c                4   AINV preconditioning
c                if iprec=3  ntprec=nterm
c                        =0-2 ntprec is not used
c imax         = maximum # of iterations
c nequ         = dimension of system matrix
c nterm        = # of nonzero elements in system matrix
c ntermp       = # of nonzero elements in preconditioning matrix
c ndir         = # of equations to exclude in the calculation of the 
c                residual
c info         = 0 normal exit (convergence achieved)
c              = 1 imax reached
c iter         = final number of iteration to convergence
c ix0          = 0 --> x0 = 0                                  
c              = 1 --> x0 = K^{-1}b                            
c iexit        = 0 exit on absolute residual               
c              = 1 exit on |r_k|/|b|                       
c
c Integer arrays
c **************
c
c noddir(ndir) = equation numbers to be excluded in the calculation of
c                the residual
c ia(nequ+1)   = pointer to diagonal element of row for system matrix
cxcx ja(nterm)    = column index for system matrix
c ja(N1*N)     = column index for system matrix
c iap(nequ+1)  = pointer to diagonal element of row for AINV matrix
c                (used only if iprec=4)
c jap(ntermp)  = column index for AINV matrix (used only if iprec=4)
c
c
c Real*8 scalars
c **************
c
c tolcg        = tolerance for relative residual
c bnorm        = norm of right hand side
c resini       = norm of initial residual (|b-Ax_0|/|b|)
c resnorm      = norm of residual at final iteration (|b-Ax_k|/|b|)
c
c
c Real*8 arrays
c *************
c 
c sysmat(nterm)= coefficient matrix
c prec(ntermp) = preconditioning matrix (triangular factor for iprec=3,4)
c rhs(nequ)    = right hand side
c sol(nequ)    = solution vector (xk)
c scr(3*nequ)  = scratch array
c res(nequ)    = PCG residual (rk)
c pk(nequ)     = PCG direction (pk)
c
c Logical scalars
c ***************
c
c exit_test    = self explanatory
c
      implicit none

      integer   i,iaxp,ipres,iainv
      
      integer   iprt,iexit,iprec,imax,nequ,nterm,ntermp,ndir
      integer   info,iter

cxcx  integer   noddir(ndir),ia(nequ+1),ja(nterm)
cxcx  integer   iap(nequ+1),jap(ntermp)
      integer   noddir(*),ia(*),ja(*)
      integer   iap(*),jap(*)

      real*8    tolcg,bnorm,resini,resnorm

      real*8    beta,alpha,zero

cxcx  real*8    sysmat(nterm),prec(ntermp)
cxcx  real*8    rhs(nequ),sol(nequ),res(nequ),pk(nequ),scr(3*nequ)
      real*8    sysmat(*),prec(*)
      real*8    rhs(*),sol(*),res(*),pk(*),scr(*)

      real*8    normres,dnrm2,ddot

      logical   exit_test

      parameter (zero=0.d0)


c
c  intialize pointers to scratch array scr
c
      iaxp=1
      ipres=iaxp+nequ
      iainv=ipres+nequ

c
c normal exit flag
c
      info = 0
      exit_test = .false.
      iter=0
      if (iprt.gt.0) write(iprt,'(a)') '      iter        resnorm'
      bnorm=normres(nequ,ndir,noddir,rhs,scr(iainv))
c
c handles case of zero RHS
c
      call zerorhs(nequ,info,ndir,bnorm,resnorm,sol,rhs,exit_test)
c
c  calculate initial residual
c
      call ressym2(nequ,nterm,ndir,ia,ja,noddir,sysmat,rhs,sol,res)
      if (iexit.eq.0) then
         resini=normres(nequ,ndir,noddir,res,scr(iainv))
      else
         resini=normres(nequ,ndir,noddir,res,scr(iainv))/bnorm
      endif

      do while (.not. exit_test)
c
         iter=iter+1
c
c  calculates PRECO.r_k+1
c     
         call precond(iprec,nequ,ntermp,iap,jap,prec,res,scr(iainv),
     1                scr(ipres))

c
c  calculates \beta_k
c
         if (iter.eq.1) then
            beta=zero
         else
            beta=-ddot(nequ,scr(ipres),1,scr(iaxp),1)/
     1            ddot(nequ,pk,1,scr(iaxp),1)
         end if
c
c  calculates p_k+1:=SCR(ipres)+beta*p_k
c
         call dxpay(nequ,scr(ipres),1,beta,pk,1)
c
c  calculates \alpha_k
c
         call axbsym(nequ,nterm,ia,ja,sysmat,pk,scr(iaxp))
         alpha = ddot(nequ,pk,1,res,1)/ddot(nequ,pk,1,scr(iaxp),1)
c
c  calculates x_k+1 and r_k+1
c
         call daxpy(nequ,alpha,pk,1,sol,1)
         call daxpy(nequ,-alpha,scr(iaxp),1,res,1)
c     
c  set up loop control
c
         if (iexit.eq.0) then
            resnorm=normres(nequ,ndir,noddir,res,scr(iainv))
         else
            resnorm=normres(nequ,ndir,noddir,res,scr(iainv))/bnorm
         endif
         if (iprt.gt.0) write(iprt,'(i10,e15.7)') iter,resnorm
         exit_test = (iter.gt.imax .or. resnorm.le.tolcg)
      end do

      if(iter.ge.imax) info=1

      return
      end
