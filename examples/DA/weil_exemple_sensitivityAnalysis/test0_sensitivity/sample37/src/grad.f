       subroutine  grad(nmax,nequ,nterm,itmax,ia,ja,coef,prec,
     1                  tnoto,tol,xnew,rmax)
C
C
       implicit none
       integer nequ,nterm,itmax,nmax
       integer i, iter
       integer ia(nequ+1), ja(nterm)

       real*8  coef(nterm), prec(nterm)
       real*8  tnoto(nequ), xold(nequ), xnew(nequ)
       real*8  rold(nequ), rnew(nequ)
       real*8  pold(nequ), pnew(nequ)
       real*8  hr(nequ),hp(nequ)
       real*8  kh(nequ),h_res(nequ), rmax

       real*8  tol,resnorm,norma,vettvett,norma1,norma2
       real*8  alpha, beta,nrm

       
       call lsolve(nequ,nterm,ia,ja,prec,tnoto,xold)
       call matvett(nmax,nequ,nterm,ia,ja,coef,xold,rold)
       do i=1,nequ
         rold(i)=tnoto(i) -rold(i)
       end do
       call lsolve(nequ,nterm,ia,ja,prec,rold,pold)
       iter=0
       resnorm=2*tol
       write(98,*) 'CONVERGENZA DEL GCM'
       do while ((iter.lt.itmax).and.(resnorm.ge.tol))
         iter = iter+1
         alpha= vettvett(nequ,pold,rold)
         call matvett(nmax,nequ,nterm,ia,ja,coef,pold,hp)
         alpha=alpha/vettvett(nequ,pold,hp)
         do i=1,nequ
              rnew(i)=rold(i) -alpha*hp(i)
         end do
         call lsolve(nequ,nterm,ia,ja,prec,hp,kh)
         beta=-vettvett(nequ,rnew,kh)/vettvett(nequ,pold,hp)
         call lsolve(nequ,nterm,ia,ja,prec,rnew,h_res)
         do i=1,nequ
            pnew(i) =h_res(i) +beta*pold(i)
            xnew(i)= xold(i) +alpha*pold(i)
         end do
         resnorm=norma2(nequ,rnew)/norma2(nequ,tnoto)*(1.e-09*rmax)
c        resnorm=norma2(nequ,rnew)/norma2(nequ,tnoto)
         write(98,101) iter, resnorm
 101     FORMAT(1x,'ITER:  ',I5,2x,'RESIDUO:  ',E16.5)
         do i=1,nequ
            pold(i) = pnew(i)
            xold(i) = xnew(i)
            rold(i) = rnew(i)
         end do
      end do
c     call matvett(nmax,nequ,nterm,ia,ja,coef,xnew,rnew)
c     do i=1,nequ
c          rnew(i)=tnoto(i)-rnew(i)
c     end do
c     resnorm=norma(nequ,rnew)/norma(nequ,tnoto)
c     write(*,*) 'residuo vero', resnorm

      return 
      end
