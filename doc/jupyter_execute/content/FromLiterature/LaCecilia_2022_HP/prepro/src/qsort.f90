!-----------------------------------------------------------------------------
!
! QSORT
!
! Quicksort
!
!-----------------------------------------------------------------------------

  subroutine qsort(n,arr,brr)

  integer(kind=4) n,M,NSTACK
  real(kind=8) arr(n)
  real(kind=4) brr(n)
  parameter (M=7,NSTACK=50)

! Sorts an array arr(1:n) into ascending order using Quicksort, while
! making corresponding rearrangement of the array brr(1:n).

  integer(kind=4) i,ir,j,jstack,k,l,istack(NSTACK)
  real(kind=8) a,b,temp
  jstack=0
  l=1
  ir=n

! Insertion sort when subarray small enough.

1 if (ir-l.lt.M) then
     do j=l+1,ir
        a=arr(j)
        b=brr(j)
        do i=j-1,1,-1
           if (arr(i).le.a) goto 2
           arr(i+1)=arr(i)
           brr(i+1)=brr(i)
        end do
        i=0
2       arr(i+1)=a
        brr(i+1)=b
     end do
     if (jstack.eq.0) return
     ir=istack(jstack) ! Pop stack and begin a new round of partiti
     l=istack(jstack-1)
     jstack=jstack-2
  else

! Choose median on left, center and right elements as partitioning eleme
! Also rearrange so that a(l+1) .leq. a(l) .leq. a(ir).

     k=(l+ir)/2
     temp=arr(k)
     arr(k)=arr(l+1)
     arr(l+1)=temp
     temp=brr(k)
     brr(k)=brr(l+1)
     brr(l+1)=temp
     if (arr(l+1).gt.arr(ir)) then
        temp=arr(l+1)
        arr(l+1)=arr(ir)
        arr(ir)=temp
        temp=brr(l+1)
        brr(l+1)=brr(ir)
        brr(ir)=temp
     end if
     if (arr(l).gt.arr(ir)) then
        temp=arr(l)
        arr(l)=arr(ir)
        arr(ir)=temp
        temp=brr(l)
        brr(l)=brr(ir)
        brr(ir)=temp
     end if
     if (arr(l+1).gt.arr(l)) then
        temp=arr(l+1)
        arr(l+1)=arr(l)
        arr(l)=temp
        temp=brr(l+1)
        brr(l+1)=brr(l)
        brr(l)=temp
     end if

! Initialize pointers for partitioning.

     i=l+1
     j=ir
     a=arr(l) ! Partitioning element.
     b=brr(l)
3    continue ! Beginning of innermost loop.
     i=i+1 ! Scan up to find element .gt. a.
     if (arr(i).lt.a) goto 3
4    continue
     j=j-1 ! Scan down to find element .lt. a.
     if (arr(j).gt.a) goto 4
     if (j.lt.i) goto 5 ! Pointers crossed. Partitioning complete.
     temp=arr(i) ! Exchange elements of both arrays.
     arr(i)=arr(j)
     arr(j)=temp
     temp=brr(i)
     brr(i)=brr(j)
     brr(j)=temp
     goto 3 ! End of the innermost loop.
5    arr(l)=arr(j) ! Insert partitioning element in both arrays.
     arr(j)=a
     brr(l)=brr(j)
     brr(j)=b
     jstack=jstack+2

! Push pointers to larger subarray on stack, process smaller subarray
! immediately.

     if (jstack.gt.NSTACK) then
        write(6,*) 'NSTACK too small!'
        stop
     end if
     if (ir-i+1.ge.j-l) then
        istack(jstack)=ir
        istack(jstack-1)=i
        ir=j-1
     else
        istack(jstack)=j-1
        istack(jstack-1)=l
        l=i
     end if
  end if
  goto 1

  end subroutine qsort
