! Chase Funkhouser
! This module provides subroutines for setting up radius and force grids and constants
module m_grids
use m_defs, only: dp
integer :: nrad
real(dp) :: rmin,rmax,cs_squared,gm,temp_const,adiabatic_index,m_planet,rho_crit
real(dp), allocatable :: rad(:),force(:)

contains


subroutine make_grids
use m_defs, only: kboltz,cgrav,mjup,mproton
implicit none
integer :: i

namelist /grid_nml/ nrad,rmin,rmax
namelist /const_nml/ temp_const,adiabatic_index,m_planet,rho_crit
open (12,file='input.nml')
read (12,nml=grid_nml)
read (12,nml=const_nml)
close(12)

allocate( rad(nrad) )

! make radius grids linearly spaced
do i=1,nrad
  rad(i) = (rmin + (rmax-rmin)*(i-1.0)/(nrad-1.0))
end do

! set up sound speed
cs_squared = adiabatic_index*temp_const*kboltz/mproton

! set GM value using m_planet in units of mjup
gm = cgrav*m_planet*mjup

end subroutine make_grids


subroutine make_force
implicit none
integer :: i
real(dp) :: a0,a1,a2,r


namelist /force_nml/ a0,a1,a2
open (12,file='input.nml')
read (12,nml=force_nml)
close(12)

allocate( force(nrad) )

force=0.d0

! set force as second degree polynomial in r
do i=1,nrad
  r=rad(i)
  force(i)=a0+a1*r+a2*r**2
enddo

end subroutine make_force


end module m_grids
