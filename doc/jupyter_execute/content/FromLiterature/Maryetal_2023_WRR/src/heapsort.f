      subroutine heapsort(face,n)
c
c  Ordina in senso crescente il vettore intero face
c      
      implicit none
      integer  node,i,j,k,ik,jk
      integer  n 
      integer  face(5,n)
      logical  conf_fac
c 
      do node = 2,n
         i = node
         j = i/2
         do while((i.ne.1).and.(conf_fac(face(1,j),face(1,i))))
            call swap_t(face(1,j),face(1,i))
            i = j
            j = i/2
         end do
      end do
c
      do i = n,2,-1
         call swap_t(face(1,i),face(1,1))
         k = i-1
         ik = 1
         jk = 2
         if ((k.ge.3).and.(conf_fac(face(1,2),face(1,3))))  jk = 3
         do while ((jk.le.k).and.(conf_fac(face(1,ik),face(1,jk))))
            call swap_t(face(1,jk),face(1,ik))
            ik = jk
            jk = ik*2
            if ((jk+1.le.k).and.(conf_fac(face(1,jk),face(1,jk+1))))
     1         jk = jk+1  
         end do
      end do
c
      return 
c
      end
