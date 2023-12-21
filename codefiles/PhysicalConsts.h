#ifndef _CONST_LIST_H
#define _CONST_LIST_H

/// List of physical constants required by the code

#define G_Newt 6.67428e-11 //usual units  - //6.2565e5 //bar/(g cm^-3) (R_E/M_E)
#define M_Earth 5.972e24 //kg
#define R_Earth 6.371e6 //m
#define Sig_Stefan 5.670374419e-8
// also recall: ac = 4*Sig_Stefan
#define boltz_R 8300
#define R_gas 8.3145

#define yrs_to_s 3.154e7

//bulk perovskite equation of state, Seagar(2007)
#define rho0 4100.0
#define aK 0.00692e9
#define Tcrit 300
#define perov_n 0.541
#define perov_C 0.00161


#endif