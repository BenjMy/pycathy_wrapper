!-----------------------------------------------------------------------------
!
! WBB
!
! This program defines the structure of the binary files basin_b/basin_i
! using the ascii DEM file <demfile>, which contains the catchment cell 
! elevations. The following reference system is considered:
!
!              ^
!              |
!              |
!    j=[1,M],y |
!              |
!              |
!             -|--------------->
!                  i=[1,N],x
!
!-----------------------------------------------------------------------------

program wbb

call wbb_sr

end program wbb

