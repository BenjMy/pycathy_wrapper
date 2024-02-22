C
C************************ LSOLVE ************************************
C
      subroutine lsolve(nequ,nterm,ia,ja,lmat,vec,pvec)
c
c  forward and backward substitution for the application
c  of the Choleski preconditioner
c
c   pvec:=(LL^T)^-1*vec (solves the system LL^T*pvec=vec)
c
      implicit none  

      integer  nequ,nterm
      integer  k,n1,mm,i,j,m
cxcx  integer  ia(nequ+1),ja(nterm)
cxcx  real*8   lmat(nterm),pvec(nequ),vec(nequ)
      integer  ia(*),ja(*)
      real*8   lmat(*),pvec(*),vec(*)
      real*8   a,zero

      parameter (zero=0.d0)

      do k = 1,nequ
         pvec(k) = zero
      end do
      do k = 1,nequ
         i = ia(k)
         j = ia(k+1) - 1
         pvec(k) = (vec(k) - pvec(k))/lmat(i)
         do m = i+1,j
            pvec(ja(m)) = pvec(ja(m)) + lmat(m)*pvec(k)
         end do
      end do
      do k = 1,nequ 
         n1 = nequ-k+1
         a = zero
         i = ia(n1)
         j = ia(n1+1) - 1
         do m = i,j-1
            mm = j-m+i
            a = a + lmat(mm)*pvec(ja(mm))
         end do
         pvec(n1) = (pvec(n1) - a)/lmat(i)
      end do

      return
      end
