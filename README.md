---
title: ASTRO 5470 Final Project
author: Chase Funkhouser
date: May 10, 2024
---

# Background

This program attempts to solve the boundary value isothermal wind problem with an additional force term. The momentum equation for a steady-state radial wind is given by:

$$ v\frac{dv}{dr}=-\frac{GM}{r^2}+\frac{1}{\rho}\frac{dP}{dr}+g_l $$

Where $v$ is the velocity at radius $r$, $M$ the planet mass, $P$ the pressure and $g_l$ the additional acceleration term. This is supplemented by the law of conservation of mass (continuity equation):

$$ \dot{M}=4\pi r^2v\rho=\text{constant} $$

We further assume that the wind is isothermal with sound speed $c_s$ that satisfies the equation of state $P=\rho c_s^2$. With an adiabatic index $\gamma$ we have $c_s=\sqrt{\gamma kT/m}$ with temperature $T$ and particle mass $m$. Then the final equation can be written in standard form as:

$$ \frac{dv}{dr}=\frac{-GM+2c_s^2+g_lr^2}{r^2(v^2-c_s^2)/v} $$

At the sonic point where $v=c_s$ the derivative is singular, so we apply the regularity condition that the numerator also be zero, giving a critical radius at $r_s$ given by:

$$ -GM + 2c_s^2r_s+g_lr_s^2=0 $$

The boundary value problem then is to solve for the transsonic solution which passes smoothly through the sonic point, starting at subsonic velocity below the critical radius and supersonic above the critical radius. In the case where $g_l=0$, this is the usual Parker solar wind problem. The addition of $g_l$ effectively shifts the critical point away from the usual Parker sonic point at $GM/2c_s^2$. Since the equation defining the velocity is nonlinear, solution of this problem requires a numerical integration.


# Methods

## Rootfinding with Bracketing

Before the code can integrate from the critical point, it has to determine where the critical point actually is. Since the critical point is at the root of a polynomial-ish expression (where $g_l$ can be any smooth function of radius this won't be the case) we use a simple rootfinding scheme. We simple move through each radius interval, determining at each step if the numerator of the derivative switches sign (we can't find an exact zero numerically, so this is close enough). The left side (smaller radius) index is chosen than as the radius of the critical point. With a large enough fidelity in radius and a simple enough and smooth $g_l$, this method should quickly find a root.

Note that the root must lie within the domain of the problem (between `rmin` and `rmax`), and that this method only identifies the first root, which in this case is the smallest root. If the form of $g_l$ is sufficiently complex that there are multiple roots, they will not be found, and the subsequent integration may then fail.


## Runge-Kutta 4-5 Integration

Since the solution to the boundary value wind problem requires a numerical integration, this code uses a Runge-Kutta 4-5 integration scheme, with adaptive stepsizing using step doubling. Effectively, for each radius interval the code computes the regular Runge-Kutta 4th order step across the full interval, as well as two steps across half the interval. The halved steps combined give a 5th order accurate solution which should be more accurate than the 4th order step. We compare the result from the 4th order step against the 5th order, and if the error is greater than the error limit `eps`, the code halves the full stepsize and tries again. Once the error is small enough, the code steps forward this stepsize. Since the stepsize could be smaller than the distance between consecutive radius bins, the code may repeat this process until the final achieved radius passes the next desired radius bin, at which point it stores the final value of the velocity it determined.

To determine the velocity over the full interval, the code starts by integrating down from the critical point to the minimum radius `rmin`. Since the rootfinding scheme always finds a critical radius below the "true" value, this radius is used to start the integration downward, with a small velocity offset down from the sound speed to avoid the singularity. Once this is finished, the code integrates outwards from the first radius bin above the critical radius (which will be above the "true" critical radius) to the maximum radius `rmax`, again avoiding the singularity by starting with a small velocity offset above the sound speed.


# Implementation

## defs.f90

This module provides constants and parameter definitions for the rest of the code.


## grids.f90

This module provides subroutines for setting up the radius and force grids. 

### make_grids
The subroutine `make_grids` sets up the radius grid by pulling parameters `rmin` and `rmax` from the input file `input.nml` under the `grid_nml` namelist. The radius grid is linear, starting at `rmin` and ending at `rmax`, with a total number of radius points `nrad` also set in the input file. The resulting grid is stored in the `rad` array. Once the radius grid is built, the code further defines important variables from the input file under than `const_nml` namelist. These include the isothermal temperature `temp_const` in Kelvin, the adiabatic index `adiabatic_index`, the planet mass `m_planet` in units of jupiter masses, and the density and the critical radius `rho_crit` in g/cm^3. Using these parameters the subroutine sets the sound speed squared `cs_squared` and the gravitational acceleration contribution `gm`. 

### make_force
The subroutine `make_force` creates the grid of additional force values in the array `force`. Currently this pulls input parameters from the `force_nml` namelist. The force is set to a second degree polynomial in radius:

$$ g_l=a_0+a_1r+a_2r^2 $$

Note that the `a1` and `a2` constants have units of force/radius and force/radius^2, respectively. These should be carefully chosen in connection with the radius domain used by the code to ensure that the force is not too large that the critical point is shifted outside of the domain of integration.


## integrate.f90

This module contains subroutines for integrating the velocity equation using a Runge-Kutta 4-5 method with adaptive stepsizing. 

### run_integrate
The subroutine `run_integrate` performs the majority of the work by implementing the adaptive stepsizing using step doubling. It takes as inputs the start and ending radius indices `istart` and `iend`, the initial velocity `v0`, and the maximum allowed integration error `eps`, and writes outputs to the array `vout`. For each radius bin `rad(i)`, the subroutine call `rk4` once over the full stepsize `dr`, and twice over half the stepsize to get 4th and 5th order accurate results, respectively. The results are compared in `err`; if `err` is less than the predefined maximum `eps` (provided in `solve.f90`) then the stepsize is halved and the process repeated. Once the error is within `eps`, the integrator actually takes a step forward the resulting stepsize. This process continues until the radius has surpassed that of the next radius bin `rad(i+incr)`, at which point it saves the velocity in `vout(i+incr)`. The maximum number of step doublings is set in the code by the `itermax` variable as 10; if this is surpassed the code will output a warning. The `incr` variable is set to either 1 or -1, depending on if the integration should proceed outwards (to larger radius) or inwards (to smaller radius).

### rk4
The subroutine `rk4` takes one 4th order Runge-Kutta step at radius `r` and stepsize `dr`. The input `g` is the force at the given radius, used in the call to `derivs` to calculate the derivative term. The initial velocity is provided in the input `v`, and the ouput velocity in `vout`.

### derivs
The subroutine `derivs` evaluates the derivative of the velocity at radius `r` and velocity `v` given force `g`. The result is stored in `dvdr`. Note that this routine can produce a singularity when the velocity is too close to the sound speed.


## solve.f90

This module provides subroutines for finding the critical radius and solving the resulting boundary value problem.

### initialize
The subroutine `initialize` sets parameters for the integration given in the `integ_nml` namelist. Currently the only parameter is `eps`, the maximum error allowed in each RK45 step. 

### find_rcrit
The subroutine `find_rcrit` finds the critical radius through rootfinding by bracketing. The equation to be solved for the critical radius is:

$$ g_lr^2+2c_s^2r-GM=0 $$

Where the array `force` stores the values of $g_l$ at each radius `rad`. The code loops through each radius interval and determines if the function flips sign within the interval by looking at the endpoints. Once a root is found, the code exits and returns the root in the variable `rcrit`. If no root is found, the code stops the program and returns an error message. The subroutine also prints a comparison of the critical radius with the usual Parker sonic radius.

### solve_bvp
The subroutine `solve_bvp` solves the boundary value problem by integrating from the critical radius `rcrit` to the boundaries `rmin` and `rmax`. It first finds the closest indexed radius to the critical radius as `icrit`, then integrates inward from `icrit` to `iend` at 1, the index for `rmin`. The initial velocity is set to the sound speed times `1.0-1.d-2`, a small deviation below 1, to ensure that the code bypasses the singularity.

A similar technique is used for the integration outwards in the second step, which integrates from `icrit+1` out to `nrad` given an initial velocity that is `1.0+1.d-2` above the sound speed. The result should be a velocity array `finvel` that is subsonic below the sonic point and supersonic above.


## dump.f90

This module provides subroutines for writing output data. The subroutine `dump_all` writes the final velocity profile from `finvel` at every radius to the output file `final_velocity.data`. It also writes all of the important setup parameters into a file `setup.data`. The subroutine `density_mdot` calculates the density profile from the velocity profile, as well as the value for $\dot{M}$ in `mdot`. These are written to the file `density_mdot.data`.


## main.f90

This module provides the main program for the code. First the grids `rad` and `force` are created using the `make_grids` and `make_force` subroutines. Then the `eps` parameter is set via the `initialize` routine. Then the solivng begins by finding the critical radius using `find_rcrit`. Once the root is found, `solve_bvp` integrates outwards to both boundaries from the critical radius. The resulting velocity profile is written to an output in the subroutine `dump_all`, and the subroutine `density_mdot` writes the density and $\dot{M}$ values out.


## Makefile

The makefile provides a quick way to compile the code and clean the directory. To compile the code, simply type `make` in the terminal within the code directory. The compiler will produce an executable file venerably named `a.out` which can be run via `./a.out`. The directory can also be cleaned using `make clean`, which will remove all object files and `a.out`.


## input.nml

This namelist provides input parameters that the code will use. The namelists are separated into four blocks.

### grid_nml
This namelist provides parameters for the radius grid:

- `nrad`: integer, number of radius bins to use. Higher numbers result in greater accuracy but longer computation time.
- `rmin`: float, minimum radius to use in the domain, in cm.
- `rmax`: float, maximum radius to use in the domain, in cm.

### const_nml
This namelist provides various constants to be used by the code:

- `temp_const`: float, isothermal temperature to use, in Kelvin.
- `adiabatic_index`: float, value of the adiabatic index.
- `m_planet`: float, planet mass in units of Jupiter masses.
- `rho_crit`: float, density at the critical radius in g/cm^3.

### integ_nml
This namelist provides parameters for the integration:

- `eps`: float, maximum allowed relative error in each Runge-Kutta integration step.

### force_nml
This namelist provides parameters for the force function:

- `a0`: float, value of the constant force contribution, as an acceleration in cgs.
- `a1`: float, value of the force contribution proportional to radius, as an acceleration/radius in cgs.
- `a2`: float, value of the force contribution proportional to radius^2, as an acceleration/radius^2 in cgs.


## plots.py

This Python script provides functionality for plotting the results from the program output. The three output files `setup.data`, `final_velocity.data`, and `density_mdot.data` are read in using the respective functions `read_setup_file`, `read_vel_file`, and `read_density_file`. The script always requires a `setup.data` file to read in the same directory. However, the data file is read from the command line argument preceded by the flag `-i`. The type of plot is also specified at the command line using the `-p` flag. The current possible plots are `velocity` to plot the velocity from `final_velocity.data`, `density` to plot density from `density_mdot.data`, and `mdot` to plot the mass loss rate from `density_mdot.data`. These are each plotted through the functions `plot_vel`, `plot_dens`, and `plot_mdot`, respectively.

The `velocity` plot will show the velocity in solid blue, as well as the sound speed in solid red, and the critical radius in dotted blue. Note that the velocity curve ought to pass through the sound speed curve at the critical radius. For analysis, the Parker radius is also provided in dashed green.

The `density` plot will show the density curve in solid blue, with a red horizontal line showing the value of the density at the critical point (the critical density), and the critical radius is shown as a dotted blue vertical line. The density curve, like the velocity curve, ought to move smoothly through the critical radius.

The `mdot` plot is a test for the accuracy of the code. Using the continuity equation, the value of `mdot` through the radius domain should be constant. If the code is workign properly, the resulting plot will show a horizonal blue line as the mass loss rate.

