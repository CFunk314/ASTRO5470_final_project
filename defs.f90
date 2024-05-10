! Chase Funkhouser
module m_defs
implicit none

integer, parameter :: longint = selected_int_kind(16)
integer, parameter :: dp = kind(1.d0)
real (dp), parameter :: pi=2.d0*asin(1.d0)
real (dp), parameter :: twopi=2.d0*pi

real (dp), parameter :: cgrav=6.67429d-8
real (dp), parameter :: hplanck=6.62606896d-27
real (dp), parameter :: clight=2.997924589d10
real (dp), parameter :: clightkms=2.997924589d5
real (dp), parameter :: kboltz=1.3806504d-16
real (dp), parameter :: charge=4.80320427d-10
real (dp), parameter :: abohr=0.52917720859d-8
real (dp), parameter :: melectron= 9.10938215d-28
real (dp), parameter :: mproton=1.672621637d-24
real (dp), parameter :: amu=1.660538782d-24
real (dp), parameter :: eV = 1.602176487d-12
real (dp), parameter :: angstrom = 1.d-8
real (dp), parameter :: kms=1.d5

real (dp), parameter :: au = 1.49597870700d13
real (dp), parameter :: msun = 1.9884d33
real (dp), parameter :: rsun = 6.9551d10
real (dp), parameter :: mjup = 1.8987d30
real (dp), parameter :: rjup = 7.1492d9

end module m_defs
