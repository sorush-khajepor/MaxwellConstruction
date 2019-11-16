#  Maxwell Construction Solver

This solver finds the saturation pressure, vapor and liquid density at a specific temperature for Van der Waals (VW), Peng-Robinson (PR), Redlich-Kwong (RK), Soave-Redlich-Kwong (SRK), Carnahan-Starling (CS). It can also produce results for a range of temperature with double precision accuracy.

prerequisite
CMake
Gfortran >= 5.0
gnuplot (optional)

 build on windows
 mkdir build
 cd build
 cmake -G "MinGW Makefiles" ..
 mingw32-make
 MaxwellConstruction.exe

 build on linux
 mkdir build
 cd build
 cmake ..
 make
 ./MaxwellConstruction
 
 Once successfully compiled, and run, you can change settings in file ./examples/main.f90
 
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
  
East peasy :)
