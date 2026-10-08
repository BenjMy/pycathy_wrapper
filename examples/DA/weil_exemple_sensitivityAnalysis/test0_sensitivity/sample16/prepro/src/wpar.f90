!-----------------------------------------------------------------------------
!
! WPAR
!
! This program reads the parameter file <parfile> and writes this file in
! a well-formatted form. Parameters of the original <parfile> should be
! provided in each line after the "=" mark. The rewritten <parfile> will
! add parameter descriptions and a formatted parameter structure as 
! specified in the module mpar.
!
!-----------------------------------------------------------------------------

  program wpar

  use mpar

  implicit none

  call rparfile('hap.in')
  call wparfile('hap.in')

  end program wpar
