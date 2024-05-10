
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
real(dp) :: a0,a1,a2,r

namelist /grid_nml/ nrad,rmin,rmax
namelist /const_nml/ temp_const,adiabatic_index,m_planet,rho_crit
namelist /force_nml/ a0,a1,a2
open (12,file='input.nml')
read (12,nml=grid_nml)
read (12,nml=const_nml)
read (12,nml=force_nml)
close(12)

allocate( rad(nrad),force(nrad) )

! make radius grids linearly spaced
do i=1,nrad
  rad(i) = (rmin + (rmax-rmin)*(i-1.0)/(nrad-1.0))
end do

! set up sound speed
cs_squared = adiabatic_index*temp_const*kboltz/mproton

! set GM value using m_planet in units of mjup
gm = cgrav*m_planet*mjup

! set up force grid as a second order polynomial in r
do i=1,nrad
  r=rad(i)
  force(i) = a0+a1*r+a2*r**2
enddo

end subroutine make_grids


end module m_grids
