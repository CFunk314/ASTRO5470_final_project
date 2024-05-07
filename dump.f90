module m_dump
use m_defs, only: dp

contains


subroutine dump_all
use m_defs, only: kms
use m_grids
use m_solve
implicit none
integer :: i
character*15 :: c_nrad

write(c_nrad,"(i0.3)") nrad

open(15,file='final_velocity_nrad'//trim(c_nrad)//'.data')
write(15,"(a15,2x,i4)") "nrad =", nrad
write(15,"(2(2x,a15))") "rad","vel"
do i=1,nrad
  write(15,"(2(2x,es15.8))") rad(i),finvel(i)/kms
end do
close(15)

end subroutine dump_all


end module m_dump
