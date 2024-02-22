c
c*************************  ORDN  ****************************
c  ordinamento crescente di un vettore di interi
      subroutine ordn(lvec,ivec)

      implicit none

      integer lvec,itemp
      integer ivec(lvec),j,k

      do j = 1,lvec-1
         do k = j+1,lvec
            if(ivec(j).gt.ivec(k)) then
               itemp = ivec(j)
               ivec(j) = ivec(k)
               ivec(k) = itemp  
            end if
         end do
      end do

      return
      end
