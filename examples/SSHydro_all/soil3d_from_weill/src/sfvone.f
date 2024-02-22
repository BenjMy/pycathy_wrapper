      subroutine sfvone(iin,iterm,iout,iprt,n,nsfmax,nnsfmx,
     1                  time,nsf,nsfnum,nsfnod,sfv,sfvnum,
     2                  sfvnod,sfvtim,ifsf,ifsfp,z)
c
      implicit none
      integer  i,j,inod1,inod2
      integer  iin,iterm,iout,iprt
      integer  n,nsfmax,nnsfmx,nsf
      integer  ifsf(*),ifsfp(*) 
      integer  nsfnum(nsfmax),nsfnod(nsfmax,nnsfmx)
      integer  sfv(2),sfvnum(2,nsfmax),sfvnod(2,nsfmax,nnsfmx)
      real*8   time,sfvtim(2)
      real*8   z(n)
c

      do i=1,n
         ifsf(i)=0
      end do 
      read(iin,*) sfvtim(1)
      read(iin,*) sfv(1)
      if (sfvtim(1).gt.time) go to 100
c
c  unit iin7 input (seepage face bc's)
c
      if (iprt.eq.2) write(iout,1400) sfv(1)
      if (iprt.eq.2) write(iterm,1400) sfv(1)
      if (sfv(1).gt.0) then
         do i=1,sfv(1)
            read(iin,*) sfvnum(2,i)
            read(iin,*) (sfvnod(2,i,j),j=1,sfvnum(2,i))
            if (iprt.eq.2) then
               write(iout,1410) i,sfvnum(2,i)
               write(iout,1420) (sfvnod(2,i,j),j=1,sfvnum(2,i))
               write(iterm,1410) i,sfvnum(2,i)
               write(iterm,1420) (sfvnod(2,i,j),j=1,sfvnum(2,i))
            end if
         end do
      end if
c

      nsf = sfv(1)
      do i=1,nsf
         nsfnum(i) = sfvnum(2,i)
         do j=1,nsfnum(i)
            nsfnod(i,j) = sfvnod(2,i,j)
         end do
         nsfnod(i,nsfnum(i)+1)=-9999
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
      do i=1,nsf
         do j=1,nsfnum(i)
            ifsf(nsfnod(i,j))=1
         end do 
      end do
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
  100 return
      end
