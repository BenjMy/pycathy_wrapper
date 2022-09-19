!-----------------------------------------------------------------------------
!
! RN
!
! The program RN determines the physiographic features of a drainage system
! using grid-based DEM data.
!
!-----------------------------------------------------------------------------

program rn

use mpar
use mbbio

implicit none

! read <parfile>

call rparfile('hap.in')

WRITE(6,*)

! open the binary files basin_b/basin_i

call load_dtm("basin_b","basin_i")

WRITE(6,*) 'csort I...'
call csort()
WRITE(6,*) '...completed'
WRITE(6,*)

WRITE(6,*) 'depit...'
call depit()
WRITE(6,*) '...completed'
WRITE(6,*)

WRITE(6,*) 'csort II...'
call csort()
WRITE(6,*) '...completed'
WRITE(6,*)

WRITE(6,*) 'cca...'
call cca()
WRITE(6,*) '...completed'
WRITE(6,*)

WRITE(6,*) 'smean...'
call smean()
WRITE(6,*) '...completed'
WRITE(6,*)

WRITE(6,*) 'dsf...'
call dsf()
WRITE(6,*) '...completed'
WRITE(6,*)

WRITE(6,*) 'hg...'
call hg()
WRITE(6,*) '...completed'
WRITE(6,*)

WRITE(6,*) 'saving the data in the basin_b/basin_i files...'
WRITE(6,*)

! save the file basin_b

call save_dtm("basin_b")

end program rn
