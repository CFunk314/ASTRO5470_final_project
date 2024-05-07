
module m_grids
use m_defs, only: dp
integer :: nrad
real(dp) :: rmin,rmax,cs_squared,gm
real(dp), allocatable :: rad(:)

contains


subroutine make_grids
use m_defs, only: kboltz,cgrav,mjup,mproton
implicit none
integer :: i
real(dp) :: temp_const,adiabatic_index,m_planet

namelist /grid_nml/ nrad,rmin,rmax
namelist /const_nml/ temp_const,adiabatic_index,m_planet
open (12,file='input.nml')
read (12,nml=grid_nml)
read (12,nml=const_nml)
close(12)

allocate( rad(nrad) )

! make radius grids evenly spaced in volume
do i=1,nrad
  rad(i) = (rmin**3 + (rmax**3-rmin**3)*(i-1.0)/nrad)**(1.0/3.0)
end do

! set up sound speed
cs_squared = adiabatic_index*temp_const*kboltz/mproton

! set GM value using m_planet in units of mjup
gm = cgrav*m_planet*mjup

end subroutine make_grids


end module m_grids
