module m_solve
use m_defs, only: dp
real(dp), allocatable :: initvel(:), finvel(:)

contains


subroutine initialize
use m_defs, only: kms
use m_grids, only: nrad
implicit none
real(dp) :: vel0

namelist /init_nml/ vel0
open (12,file='input.nml')
read (12,nml=init_nml)
close(12)

allocate( initvel(nrad), finvel(nrad) )

initvel=0.d0
initvel(1)=vel0*kms

end subroutine initialize


subroutine solve_bvp
use m_integrate
implicit none

call run_integrate(initvel,finvel)

end subroutine solve_bvp

end module m_solve
