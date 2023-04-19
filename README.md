#  Maxwell Construction Solver

This solver finds the saturation pressure, vapor and liquid density at a saturation temperature (Tsat) for Van der Waals (VW), Peng-Robinson (PR), Redlich-Kwong (RK), Soave-Redlich-Kwong (SRK), Carnahan-Starling (CS). It can also produce results for a range of temperature with double precision accuracy.

The code finds saturation pressure (Psat), corresponding to Tsat,  that satisfies Maxwell's construction. The algorithm in the nutshell is:
1) Assume a good Psat.
2) Solve P-Psat=0 to find vapor and liquid densities. For a saturation state, we should get three densities. The middle one is not real. 
3) Calculate Maxwell integration over the interval  [vapor density, liquid density]. This must be zero.
4) Go to one.

I created and used this library during my PhD on lattice Boltzmann method. I hope it helps others.
For more physics background see https://en.wikipedia.org/wiki/Maxwell_construction

## Is there a faster way to Find saturation state of a fluid?
Yes, instead of Maxwell construction condition, we could use the condition that the chemical potential of the liquid phase and the vapor phase should be equal. For simplicity, we can use the equality of fugacity coefficients for both phases:

φ_Liquid = φ_Vapor

This is not implemented in this code.

## Prerequisites
* CMake
* Gfortran >= 5.0
* gnuplot (optional)

 ## Build on windows
 ```bash
 mkdir build
 cd build
 cmake -G "MinGW Makefiles" ..
 mingw32-make
 MaxwellConstruction.exe
```

 ## Build on linux
 ```bash
 mkdir build
 cd build
 cmake ..
 make
 ./MaxwellConstruction
 ```
 Once successfully compiled, and run, you can change settings in file ./examples/main.f90
 
 ```fortran
!   Start reduced temperature
    Tr_st = 0.3_wp
!   Reduced temperature step size    
    dTr = 0.01_wp
!   Number of temperature steps
    Tr_steps = 80
!   Draw p-rho plot at this temperature (needs gnuplot installed)     
    Tr_plot = 0.7_wp
!   Equation of state
    allocate(srk_t::eos)
  ```
  
  Then, make and run again to get new results in a CSV file format
  
  ```bash
  make
 ./MaxwellConstruction
 ```
  
East peasy :)

## Cite this

If this code is useful to your research, please cite and read papers which I developed this code for:

```
Multipseudopotential interaction: A solution for thermodynamic inconsistency in pseudopotential lattice Boltzmann models
S Khajepor, J Wen, B Chen
Physical Review E 91 (2), 023301, year 2015

Multipseudopotential interaction: A consistent study of cubic equations of state in lattice Boltzmann models
S Khajepor, B Chen
Physical Review E 93 (1), 013303, year 2016
```



