! Chase Funkhouser
! This module provides subroutines for finding the critical radius and solving the boundary value problem
module m_solve
use m_defs, only: dp
real(dp) :: eps,rcrit
real(dp), allocatable :: finvel(:)

contains


! set the eps precision for integration from input
subroutine initialize
implicit none

namelist /integ_nml/ eps
open (12,file='input.nml')
read (12,nml=integ_nml)
close(12)

end subroutine initialize


! find the critical radius associated with the provided force, sound speed, and GM
subroutine find_rcrit
use m_grids
implicit none
integer :: i,icrit
real(dp), allocatable :: func(:)
real(dp) :: r,r_parker

allocate( func(nrad) )

! set up function for rootfinding
do i=1,nrad
  r=rad(i)
  func(i)=force(i)*r**2+2*cs_squared*r-gm
enddo

! loop through intervals until root is bracketed
do i=1,nrad-1
  if (func(i)*func(i+1) < 0) then
    print *,"root found near i,r=",i,rad(i)
    icrit=i
    exit
  endif

  if (i.eq.nrad-1) then
    print *,"no root found within domain! stopping..."
    stop
  endif
enddo

rcrit=rad(icrit)

! compare to the parker radius (no force)
r_parker=gm/2/cs_squared
print *,"r_parker=",r_parker
print *,"difference=",rcrit-r_parker

end subroutine find_rcrit


! solve the boundary value problem by integrating from rcrit outwards to both boundaries
subroutine solve_bvp
use m_grids
use m_integrate
implicit none
integer :: istart,iend,icrit
real(dp) :: v0

if ( .not. allocated(finvel) ) allocate( finvel(nrad) )
finvel=0.d0

! find the closest index to the critical radius
icrit=1+floor((nrad-1.0)*(rcrit-rmin)/(rmax-rmin))

! integrate from the critical radius down to rmin
istart=icrit
iend=1
v0=sqrt(cs_squared)*(1.0-1.d-2)
print *
call run_integrate(istart,iend,v0,finvel,eps)

! integrate from the critical radius up to rmax
istart=icrit+1
iend=nrad
v0=sqrt(cs_squared)*(1.0+1.d-2)
print *
call run_integrate(istart,iend,v0,finvel,eps)

end subroutine solve_bvp


end module m_solve
