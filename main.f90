program main
use m_defs
use m_grids
use m_solve
use m_dump
implicit none

print *,"setting up grids and initializing..."
call make_grids
call initialize

print *,"solving..."
call solve_bvp

print *,"finished solve. writing data..."
call dump_all

end program main
