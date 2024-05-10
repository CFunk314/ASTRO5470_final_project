module m_solve
use m_defs, only: dp
real(dp) :: vel0,eps,rcrit
real(dp), allocatable :: finvel(:)

contains


subroutine initialize
use m_defs, only: kms
use m_grids, only: nrad
implicit none

namelist /init_nml/ vel0,eps
open (12,file='input.nml')
read (12,nml=init_nml)
close(12)

vel0=vel0*kms

end subroutine initialize


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
    print *,"no root found! stopping..."
    stop
  endif
enddo

rcrit=rad(icrit)

r_parker=gm/2/cs_squared
print *,"r_parker=",r_parker
print *,"difference=",rcrit-r_parker

end subroutine find_rcrit


! solve the problem by integrating from rcrit outwards to both boundaries
subroutine solve_bvp
use m_grids
use m_integrate
implicit none
integer :: istart,iend,icrit
real(dp) :: v0

if ( .not. allocated(finvel) ) allocate( finvel(nrad) )
finvel=0.d0

icrit=1+floor((nrad-1.0)*(rcrit-rmin)/(rmax-rmin))

istart=icrit
iend=1
!v0=vel0
v0=sqrt(cs_squared)*(1.0-1.d-2)
print *
call run_integrate(istart,iend,v0,finvel,eps)

istart=icrit+1
iend=nrad
v0=sqrt(cs_squared)*(1.0+1.d-2)
print *
call run_integrate(istart,iend,vel0,finvel,eps)

end subroutine solve_bvp


end module m_solve
