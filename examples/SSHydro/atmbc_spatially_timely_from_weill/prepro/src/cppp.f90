!-----------------------------------------------------------------------------
!
! CPPP
!
! The Cathy PreProcessing Program includes the subroutine WBB, WPAR, CBB, RN,
! MRBB and BB2SHP. Starting from DEM data and hap.in parameters file a 
! complete set of inputs files, describing physiographic features of a drainage
! system, for CATHY simulations and GIS visualization are producted.
!
!-----------------------------------------------------------------------------

program cppp

use mpar
use mbbio

implicit none

call rparfile('hap.in')

write(6,*) 'wbb...'
call wbb_sr
write(6,*) '...wbb completed'
write(6,*)

call wparfile('hap.in')

WRITE(6,*) 'rn...'

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

write(6,*) '...rn completed'
write(6,*)

write(6,*) 'mrbb...'
call mrbb_sr
write(6,*) '...mrbb completed'
write(6,*)

write(6,*) 'bb2shp...'
call bb2shp_sr
write(6,*) '...bb2shp completed'

end program cppp
