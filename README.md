#  Maxwell Construction Solver

This solver finds the saturation pressure, vapor and liquid density at a specific temperature for Van der Waals (VW), Peng-Robinson (PR), Redlich-Kwong (RK), Soave-Redlich-Kwong (SRK), Carnahan-Starling (CS). It can also produce results for a range of temperature with double precision accuracy.

I created and used this library during my PhD on lattice Boltzmann method. I hope it helps others.
For more physics background see https://en.wikipedia.org/wiki/Maxwell_construction

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



