module m_integrate
use m_defs, only: dp

contains


! integrate using RK4 from rstart to rend given initial velocity at rstart of v0, and write output to vout
! radius is converted to units of parker radius GM/2c_s^2, and velocity to sound speed units c_s
subroutine run_integrate(istart,iend,v0,vout,eps)
use m_grids, only: rad,nrad,rmin,rmax,force
implicit none
integer :: i,incr,iter,itermax
integer, intent(in) :: istart,iend
real(dp), intent(in) :: v0,eps
real(dp), intent(out) :: vout(nrad)
real(dp) :: r,dr,v,vnext,verr,err,gi

! maximum step doublings
itermax=10

! set the increment to integrate forwards or backwards
incr=1
if (iend < istart) then
  incr=-1
endif

print *,"starting integration..."
print *,"rstart,istart=",rad(istart),istart
print *,"rend,iend=",rad(iend),iend

! set the initial velocity at rstart
v=v0
vout(istart)=v
print *,"initial velocity: ",vout(istart)

! loop RK4 with adaptive step-doubling over all indexes from istart to iend
do i=istart,iend-incr,incr
  r = rad(i)
  gi = force(i)

  ! keep looping until we pass the next radius index to save the value of v
  do while ((r-rad(i+incr))*incr < 0.d0)
    
    ! set dr; note that this resets dr after rk4 exits to prevent intervals halving continuously
    dr = rad(i+incr)-rad(i)

    do iter=1,itermax

      ! call an rk4 step over the full stepsize
      call rk4(r,v,gi,vnext,dr)
      !print *,"i,r,dr,v,vnext=",i,r,dr,v,vnext

      ! call two rk4 steps over half the stepsize to determine error
      call rk4(r,v,gi,verr,dr/2.0)
      call rk4(r+dr/2.0,verr,gi,verr,dr/2.0)
    
      ! calculate error and determine exit condition
      err = abs((vnext-verr)/vnext)
      !print *,"err=",err
      if (err < eps) then
        v=vnext
        r=r+dr
        exit
      endif

      ! error is not met, so halve step size and try again
      dr=dr/2.0

      if (iter.eq.10) then
        !print *,"itermax reached at i=",i
      endif

    enddo
    
    r=r+dr
    v=vnext

  enddo

  ! set the calculated velocity value at the indexed radius
  vout(i+incr)=v

enddo

print *,"integration finished"
return

end subroutine run_integrate


! runge-kutta 4th order, takes input r, v and step dr and stores output vout
subroutine rk4(r,v,g,vout,dr)
implicit none
integer n
real(dp) :: k1,k2,k3,k4
real(dp), intent(in) :: r,v,g,dr
real(dp), intent(out) :: vout

call derivs(r,v,g,k1)
call derivs(r,v+k1*0.5*dr,g,k2)
call derivs(r,v+k2*0.5*dr,g,k3)
call derivs(r,v+k3*dr,g,k4)

vout = v+dr*(k1+2.0*(k2+k3)+k4)/6.0

return

end subroutine rk4


! evaluates derivative of v at r; stores in dvdr
subroutine derivs(r,v,g,dvdr)
use m_grids
implicit none
real(dp), intent(in) :: r,v,g
real(dp), intent(out) :: dvdr
real(dp) :: numer,denom

numer = -gm+2.0*cs_squared*r+g*r**2
denom = r**2*(v**2-cs_squared)/v
dvdr = numer/denom

return

end subroutine derivs


end module m_integrate
