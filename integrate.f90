module m_integrate
use m_defs, only: dp

contains


subroutine run_integrate(v0,vout)
use m_grids, only: rad,nrad
implicit none
integer :: i
real(dp), intent(in) :: v0(nrad)
real(dp), intent(out) :: vout(nrad)
real(dp) :: dr,v

print *,"starting integration..."
v=v0(1)
vout(1)=v

do i=1,nrad-1
  dr = rad(i+1)-rad(i)
  call rk4(rad(i),vout(i),v,dr)
  vout(i+1)=v
end do

return

end subroutine run_integrate


subroutine rk4(r,v,vout,dr)
implicit none
integer n
real(dp) :: k1,k2,k3,k4
real(dp), intent(in) :: r,v,dr
real(dp), intent(out) :: vout

call derivs(r,v,k1)
call derivs(r,v+k1*0.5*dr,k2)
call derivs(r,v+k2*0.5*dr,k3)
call derivs(r,v+k3*dr,k4)

vout = v+dr*(k1+2.0*(k2+k3)+k4)/6.0

return

end subroutine rk4


subroutine derivs(r,v,dvdr)
use m_grids
implicit none
real(dp), intent(in) :: r,v
real(dp), intent(out) :: dvdr
real(dp) :: numer,denom

numer = -gm+2.0*cs_squared*r
denom = r**2*(v**2-cs_squared)/v
dvdr = numer/denom

return

end subroutine derivs


end module m_integrate
