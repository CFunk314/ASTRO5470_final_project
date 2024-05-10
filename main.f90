program main
use m_defs
use m_grids
use m_solve
use m_dump
implicit none

print *,"setting up grids..."
print *
call make_grids
call initialize

print *
print *,"solving..."
print *
call find_rcrit
call solve_bvp

print *
print *,"finished solve. writing data..."
print *
call dump_all
call density_mdot

print *
print *,"finished."

end program main 
