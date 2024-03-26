# README

# STEREO-D
### Structural & Thermal Evolution of Rocky Exoplanets in One-Dimension

## INSTALLATION

The program uses the **BOOST** (https://www.boost.org/) and **EIGEN** (https://eigen.tuxfamily.org) header libraries. Both are easy to install as they are simply header files. The user can either place them in a standard 'include' directory or specify their locations with the variable PDIR in 'makefile'. Either make a 'packages' directory, or change it to your own.

In addition, the line:

`cardinal_cubic_hermite() = default;` and `pchip() = default;`

must be added as a public definition in the classes `cardinal_cubic_hermite` (boost/math/interpolators/cubic_hermite.hpp) and `pchip` (boost/math/interpolators/pchip.hpp) in order for the program's implementation of the class to work.

## OPERATION
To compile all executables run 

`make all`

or to compile individual executables run `make` followed by the executable name.

To run executables I have created `bash` scripts which are in the formate `run_???.sh` which both create folders to store results in and run the executables with the required input parameters. Either modify these, create new ones or use the scripts inside them. The required input parameters for the executbales (e.g., planet mass, grid resolution) are written in the `run_???.sh` scripts and the corresponding `.cc` files in `frontcodefiles`. I also include a `submit_???.pbs` file for use on HPC.

Units are SI, unless stated otherwise.

### Executable options
`MixStruct` - with "Mix" EoS and analyitc BCs (likely the desired method)

`SeaStruct` - with Seager analyitc EoS and analyitc BCs

`StixStruct` - with "Stixrude" EoS

`AbeStruct` - with "Abe" EoS

`BC_Struct` - with "Mix" EoS and tabulated EoS

`GenStruct` - Static strucuture with "Ideal" gas EoS

`GenStruct_Tab` - `GenStruct` - but with tabulated EoS

`TimeStruct` - time evolving structure with "Ideal" gas EoS

`IrrEnvStruct2` - with "Ideal" gas EoS, analyitc EoS and only models an atmospheric envelope (has an unresolved core)

`IrrEnvStruct2_Tab` - `IrrEnvStruct` but with tabulated EoS

#### Tests

`Struct_test_boost` - Tests that can reproduce the analytically derviable strucuture of a static n=1 polytrope

`Struct_test_boost_Tab` - `Struct_test_boost` but with tabulated EoS

`Poly_Time_test_boost` - Tests that can reproduce an analytically derivable time evolution of a polytrope with special boundary conditions (see my thesis)

`Poly_Time_test_boost_Tab` - `Poly_Time_test_boost` but with tabulated EoS

`Struct_test_boost_Irr` - Tests that can reproduce the 

`Struct_test_boost_Irr_Tab` - `Struct_test_boost_Irr` but with tabulated EoS


### Equation of state options
Equations of state are included at compile time and read in at the start. The equation of state options are 
"Mix" - Selected rocky equations of state with melting included and an optional iron core (EoSs are tabulated)
"Ideal" - analytic ideal gas EoS
"Ideal_Tab" - tabulated ideal gas EoS
"Abe" - Adaptation of "Mix" where Cp is set constant. Used to try to match results from Abe's 90s papers (not worth using)
"Seager" - Simple analyitic rock equation of state from Seager paper
"Stixrude" - EoS from Stxrude for enstatite. No melting or core.

### NABLA type codes

'i' - isothermal  
'a' - adiabatic  
'c' - mixing length convection (with latent heats)  
'B' - Bower formulation (latent heats in mixing length AND F_mix)  
'g' - 'c' but with gravitational settling

### MIXING LENGTH PRESCRIPTIONS 

'c' - constant R/4  
'n' - distance to nearest boundary 

### TEMPERATURE BOUNDARY CONDITIONS

'f' - fixed  
'b' - atmosphereless black body  
'i' - Irradiated atmosphere, Guillot (2010)  
'm' - modified black body with atmpshere of optical depth \tau  
'l' - analytic fit to our grid of boundary models

### PRESSURE BOUNDARY CONDITIONS

'f' - fixed  
'o' - atmosphere with specified opacity at particular \tau surface  
'c' - P propto R^-4 to keep self-similarity of polytropes  

### Outut files

"scale" FILE CONTAINS

`M_planet P_c R_edge T_c L_edge M_core R_core`

"results" FILE CONTAINS

`m P r T L (<-scaled to max value, unscaled ->) rho nabla nabla-nabla_Ad phi T_sol T_liq Errors(on P,r,T,L)`

"fluxes" FILE CONTAINS

`L_tot L_cond L_conv L_grav L_mix`

"other" FILE CONTAINS

`Cp del nu lmix DS DV u_conv`

"bulk" FILE depends on the exact planet specifications. 
FOR "MIX" EOS IT CONTANS

`time (yrs), P_c, R_edge, T_c, L_edge, P_surf, T_surf, r_RCB, r_solidus (upper), r_liquidus, r_{critical melt frac}, r_solidus (deeper), r_{5% melt frac} (upper), r_{5% melt frac} (deeper), E_grav, E_thermal, E_total, M_{edge of grid}, M_total, Mass loss, z_edge, M_CMB, R_CMB, L_CMB, T_CMB, m_solidus (upper), m_solidus (deeper), m_liquidus, m_{critical melt frac}, P_solidus (upper), P_solidus (deeper), P_liquidus, P_{critical melt frac}, phi_outer, T_1Gpa, dt (yrs), timestep number, grid size`

"mass" FILE is a reduced version just with properties related to mass evolution. 
IT CONTAINS

`time (yrs), M_total, M_dot (lost from planet), M_dot (from day to nightside), total mass lost from dayside, R_total`

## Plotting results
In the folder `plotting` I have placed some example scripts for plotting results. These will require modification to the user's particular desires.