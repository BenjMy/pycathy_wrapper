!-----------------------------------------------------------------------------
!
! FACET
!
! The subroutine FACET calculates the aspect and the slope (positive
! downward) of the steepest direction within a given triangle (facet).
!
!-----------------------------------------------------------------------------

subroutine facet(e0,e1,e2,r,s_max_facet)

use mpar

implicit none

real(kind=REP) pi
real(kind=REP) s1,s2,e0,e1,e2,r,s_max_facet,sp,sd,s1_10,s1_11

pi=4.0_REP*datan(1.0_REP)

s1=(e0-e1)/delta_x ! delta_x=delta_y
s2=(e1-e2)/delta_x ! delta_x=delta_y
if (abs(s1).lt.epsilon(s1)) then
   if (s2.ge.0.0_REP) then
      r=+pi/2.0_REP ! valore di comodo per s1=s2=0
   else ! if (s2.lt.0.0E0) then
      r=-pi/2.0_REP
   end if
else
   r=datan(s2/s1)
end if
sp=dsqrt(s1**2.0_REP+s2**2.0_REP)
sd=(e0-e2)/dsqrt(delta_x**2.0_REP+delta_y**2.0_REP)

if (r.ge.0.0_REP.and.r.le.pi/4.0_REP.and.s1.ge.0.0_REP) then
   s_max_facet=sp
else
   if (s1.gt.sd) then
      s_max_facet=s1
      r=0.0_REP
   else
      s_max_facet=sd
      r=pi/4.0_REP
   end if
end if

return
end subroutine facet
