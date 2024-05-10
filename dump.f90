! Chase Funkhouser
! Subroutines for writing data to .data files
module m_dump
use m_defs, only: dp

contains


subroutine dump_all
use m_defs, only: kms
use m_grids
use m_solve
implicit none
integer :: i
real(dp) :: sound_speed

open(15,file='final_velocity.data')
write(15,"(a15,2x,i6)") "nrad=", nrad
write(15,"(3(2x,a15))") "rad","force","vel"
do i=1,nrad
  write(15,"(3(2x,es15.8))") rad(i),force(i),finvel(i)/kms
end do
close(15)

sound_speed = sqrt(cs_squared)/kms

open(16,file='setup.data')
write(16,*) "nrad=",nrad
write(16,*) "rmin(cm)=",rmin
write(16,*) "rmax(cm)=",rmax
write(16,*) "temp(K)=",temp_const
write(16,*) "adiabatic_index=",adiabatic_index
write(16,*) "m_planet(m_jup)=",m_planet
write(16,*) "sound_speed(kms)=",sound_speed
write(16,*) "GM(cgs)=",gm
write(16,*) "rcrit(cm)=",rcrit
close(16)

end subroutine dump_all


subroutine density_mdot
use m_defs, only: pi
use m_grids
use m_solve
implicit none
integer :: i
real(dp) :: density(nrad),mdot(nrad)

do i=1,nrad
  density(i)=rho_crit*sqrt(cs_squared)*rcrit**2/finvel(i)/rad(i)**2
  mdot(i)=4*pi*rad(i)**2*finvel(i)*density(i)
enddo

open(16,file='density_mdot.data')
write(16,"(a15,2x,i6)") "nrad=",nrad
write(16,"(3(2x,a15))") "rad","density","mdot"
do i=1,nrad
  write(16,"(3(2x,es15.8))") rad(i),density(i),mdot(i)
enddo
close(16)

end subroutine density_mdot


end module m_dump
