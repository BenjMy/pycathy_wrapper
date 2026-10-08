C
C************************* CGSOLV ***********************************
C
      subroutine cgsolv(nequ,nterm,ndir,niaux,nraux,iparm,noddir,
     1                  ia,ja,iaux,tolcg,aux,sysmat,rhs,sol)
c
c
c  driver for symmetric conjugate gradient method
c  (cfr. Gambolati 1994, page 148)
c
c  uses the symmetric CSR storage scheme for the sparse matrix
c
c  Exchange parameters
c  *******************
c
c  Integer scalars
c  ***************
c
c nequ         = dimension of system matrix
c nterm        = # of nonzero elements in system matrix
c ndir         = # of equations to exclude in the calculation of the
c                residual
c niaux        = dimension of iaux >= ntermp+nequ+1
c nraux        = dimension of aux >= ntermp+5*nequ
c                N.B. ntermp >= nterm is used when the preconditioner
c                has a number of nonzero elements larger than
c                the system matrix (e.g. AINV preconditioner)
c                Usually ntermp=nterm.
c
c Integer arrays
c **************
c
c iparm(10)    = array of control parameters (see local integer scalars)
c noddir(ndir) = equation numbers to be excluded in the calculation of
c                the residual
c ia(nequ+1)   = pointer to diagonal element of row for system matrix
cxcx ja(nterm)    = column index for system matrix
c ja(N1*N)     = column index for system matrix
c iaux(niaux)  = scratch array
c
c
c Real*8 scalars
c **************
c
c tolcg        = tolerance for relative residual
c
c
c Real*8 arrays
c *************
c
c aux(nraux)   = scratch array
c sysmat(nterm)= coefficient matrix
c rhs(nequ)    = right hand side
c sol(nequ)    = solution vector (xk)
c
c
c Local variables
c ***************
c
c Integer scalars
c ***************
c
c ntermp       = # of nonzero elements in preconditioning matrix
c iout         = output of # iterations, and different residuals, etc.
c                (<-iparm(1))
c iprt         = 0 does not print the convergence profile
c                  (n. of iteration, PCG relative residual norm)
c              > 0 prints the convergence profile in unit "iprt"
c                (<-iparm(2))
c iprec        = 0   no preconditioning
c                1   D^-1 preconditioning
c                2   D^-1-norm preconditioning
c                3   Choleski preconditioning
c                    (for compatibility with the nonsymmetric solver
c                     krysol we accept also 30,31,32)
c                4   AINV preconditioning
c                if iprec=3  ntprec=nterm
c                        =0-2 ntprec is not used
c                (<-iparm(3))
c         NB --> if iprec<0 then reuse previously calculated 
c                preconditioner
c info         = 0 normal exit (convergence achieved)
c              = 1 imax reached
c                (<-iparm(4))
c iter         = final number of iteration to convergence
c                (<-iparm(5))
c imax         = maximum # of iterations
c                (<-iparm(6))
c isol         = 0 --> SOL      is the initial solution for pcg
c              = 1 --> PREC*RHS is the initial solution for pcg
c                (<-iparm(7))
c iexit        = 0 --> exit on absolute residual                
c              = 1 --> exit on |r_k|/|b|
c                (<-iparm(10))
c
c Real*8 scalars
c **************
c
c bnorm        = norm of right hand side
c resini       = norm of initial residual (|b-Ax_0|/|b|)
c resnorm      = norm of residual at final iteration (|r_k|/|b|)
c resreal      = norm of real residual at final iteration (|b-Ax_k|/|b|)
c
      implicit none

      integer  nequ,nterm,ndir,niaux,nraux
cxcx  integer  iparm(10),noddir(ndir),ia(nequ+1),ja(nterm),iaux(niaux)
      integer  iparm(10),noddir(*),ia(*),ja(*),iaux(*)

      integer  ntermp,iout,iprt,iexit,iprec,info,iter,imax,isol
      integer  indiap,indjap,ires,ipk,iscr1,indprec

      integer  i,j

      real*8   tolcg

cxcx  real*8   aux(nraux),sysmat(nterm),rhs(nequ),sol(nequ)
      real*8   aux(*),sysmat(*),rhs(*),sol(*)
      real*8   bnorm,resini,resnorm,resreal

      real*8   normres




      iout = iparm(1)
      iprt = iparm(2)
      iprec = iparm(3)
      if(iprec.ge.30) iprec=iprec/10
      imax = iparm(6)
      isol = iparm(7)
      iexit = iparm(10)
c
c  set pointers to integer scratch array
c
      indiap = 1
      indjap = nequ + 2
c
c  set pointers to real scratch array
c
      ires = 1
      ipk = ires + nequ
      iscr1 = ipk + nequ
      indprec = iscr1 + 3*nequ
c
c  evaluates preconditioner
c
      if (iprec.ge.0)  then
         if(iprec.eq.1) then
            do i=1,nequ
               j=indprec+i-1
               aux(j)=1.0d0/sysmat(ia(i))
            end do
            ntermp=nequ
         else if (iprec.eq.3) then
            call kersh(iout,nequ,nterm,ia,ja,sysmat,aux(indprec))
            do i=1,nequ+1
               iaux(i) = ia(i)
            end do
            do i=1,nterm
               iaux(indjap+i-1) = ja(i)
            end do
            ntermp = nterm
         else 
            write(iout,*) ' Wrong iprec value, iprec=',iprec
            stop
         end if
      else
         iprec=-iprec
      end if
      if (isol.ne.0) then
         call precond(iprec,nequ,ntermp,iaux(indiap),iaux(indjap),
     1                aux(indprec),rhs,aux(1),sol)
      end if

c
c  calls the pcg routine
c
      call pcg(iprt,iexit,iprec,imax,nequ,nterm,ntermp,ndir,info,
     1         iter,noddir,ia,ja,iaux(indiap),iaux(indjap),
     2         tolcg,bnorm,resini,resnorm,sysmat,aux(indprec),
     3         rhs,sol,aux(iscr1),aux(ires),aux(ipk))

      iparm(4) = info
      iparm(5) = iter
c
c  calculates real residual
c
      call ressym2(nequ,nterm,ndir,ia,ja,noddir,sysmat,rhs,sol,
     1     aux(ires))
      if (iexit.eq.0) then
         resreal=normres(nequ,ndir,noddir,aux(ires),aux(ipk))
      else
         resreal=normres(nequ,ndir,noddir,aux(ires),aux(ipk))/bnorm
      end if

c
c  prints convergence detailes
c
      write(iout,'(2i5,3(1pe12.5),a)')
     1 info,iter,resini,resnorm,resreal,' <<symmetric solver>>'

      return
      end
