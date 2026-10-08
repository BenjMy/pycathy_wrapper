      subroutine sfvbak(nsfmax,nnsfmx,nsf,nsfnum,nsfnod,
     1                  sfv,sfvnum,sfvnod)
c
      implicit none
      integer  i,j
      integer  nsfmax,nnsfmx,nsf
      integer  nsfnum(nsfmax),nsfnod(nsfmax,nnsfmx)
      integer  sfv(2),sfvnum(2,nsfmax),sfvnod(2,nsfmax,nnsfmx)
c
      nsf = sfv(1)
      do i=1,nsf
         nsfnum(i) = sfvnum(1,i)
         do j=1,nsfnum(i)
            nsfnod(i,j) = sfvnod(1,i,j)
         end do
      end do
c
      return
      end
