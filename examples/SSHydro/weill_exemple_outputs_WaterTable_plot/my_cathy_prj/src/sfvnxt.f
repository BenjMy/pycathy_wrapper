      subroutine sfvnxt(iin,iterm,iout,iprt,n,nsfmax,nnsfmx,
     1                  nsf,nsfnum,nsfnod,sfv,sfvnum,sfvnod,
     2                  sfvtim,z)
c
      implicit none
      integer  i,j,inod1,inod2
      integer  iin,iterm,iout,iprt
      integer  n,nsfmax,nnsfmx,nsf
      integer  nsfnum(nsfmax),nsfnod(nsfmax,nnsfmx)
      integer  sfv(2),sfvnum(2,nsfmax),sfvnod(2,nsfmax,nnsfmx)
      real*8   sfvtim(2)
      real*8   z(n)
c
      do i=1,sfv(1)
         sfvnum(1,i) = sfvnum(2,i)
         do j=1,sfvnum(1,i)
            sfvnod(1,i,j) = sfvnod(2,i,j)
         end do
      end do
c
      if (iprt.ge.1) write(iout,1400) sfv(2)
      if (iprt.ge.1) write(iterm,1400) sfv(2)
      if (sfv(2).gt.0) then
         do i=1,sfv(2)
            read(iin,*) sfvnum(2,i)
            read(iin,*) (sfvnod(2,i,j),j=1,sfvnum(2,i))
            if (iprt.ge.1) then
               write(iout,1410) i,sfvnum(2,i)
               write(iout,1420) (sfvnod(2,i,j),j=1,sfvnum(2,i))
               write(iterm,1410) i,sfvnum(2,i)
               write(iterm,1420) (sfvnod(2,i,j),j=1,sfvnum(2,i))
            end if
         end do
      end if
c
      nsf = sfv(2)
      do i=1,nsf
         nsfnum(i) = sfvnum(2,i)
         do j=1,nsfnum(i)
            nsfnod(i,j) = sfvnod(2,i,j)
         end do
      end do
c
      do i=1,nsf
         do j=1,nsfnum(i)-1
            inod1 = nsfnod(i,j)
            inod2 = nsfnod(i,j+1)
            if (z(inod1).lt.z(inod2)) then
               write(iout,1500) i,j,inod1,j+1,inod2
               write(iterm,1500) i,j,inod1,j+1,inod2
               call closio
               stop
            end if
         end do
      end do
c
      sfvtim(1) = sfvtim(2)
      sfv(1) = sfv(2)
c
      read(iin,*) sfvtim(2)
      read(iin,*) sfv(2)
c
 1400 format(/,5x,'nsf  (# of seepage faces)               = ',i6)
 1410 format(  5x,'number of nodes on seepage face ',i6,'  = ',i6)
 1420 format(  5x,'node #''s : ',10i6)
 1500 format(//,' input error : elevation values not in descending',
     1          ' order on seepage face ',i6,
     2       /,4x,'nodes',i4,' (node #',i6,') and',i4,' (node #',i6,')')
c
      return
      end
