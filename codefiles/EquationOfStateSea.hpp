#ifndef _SEAGER_EQN_OF_STATE_H
#define _SEAGER_EQN_OF_STATE_H

#include <iostream>
#include <Eigen/Dense>
#include "PhysicalConsts.h"

using namespace Eigen;

void EoSrhosea(const ArrayXd P,  const ArrayXd T, ArrayXd &rho , double f)
{   
    ArrayXd T_T = (T > Tcrit).select(f*aK*(T-Tcrit) , 0);
    ArrayXd extrarho = (P >= T_T).select(pow(P - T_T , perov_n) , -pow(T_T - P  , perov_n));
    rho =  rho0 + perov_C * extrarho;  
}

void DiffrhoPsea(const ArrayXd P, const ArrayXd T, ArrayXd &diffrho , double f )
{   
    ArrayXd T_T = (T > Tcrit).select(f*aK*(T-Tcrit) , 0);
    ArrayXd bracket = (P >= T_T).select(pow(P - T_T , perov_n-1) , pow(T_T - P  , perov_n-1));
    diffrho = perov_n * perov_C * bracket ;
}

void DiffrhoTsea(const ArrayXd P, const ArrayXd T , ArrayXd &diffrho , double f )
{   ArrayXd T_T = (T > Tcrit).select(f*aK*(T-Tcrit) , 0);
    ArrayXd bracket = (P >= T_T).select(pow(P - T_T , perov_n-1) , pow(T_T - P  , perov_n-1));
    diffrho = (T > Tcrit).select(perov_n * perov_C * bracket * -f*aK , 0);
}

void Diff2rhoTTsea(const ArrayXd P, const ArrayXd T , ArrayXd &diffrho , double f)
{   ArrayXd T_T = (T > Tcrit).select(f*aK*(T-Tcrit) , 0);
    ArrayXd bracket = (P >= T_T).select(pow(P - T_T , perov_n-2) , -pow(T_T - P  , perov_n-2));
    diffrho = (T > Tcrit).select(perov_n * (perov_n -1)* perov_C * bracket * -f*aK * -f*aK, 0);
}
void Diff2rhoPTsea(const ArrayXd P, const ArrayXd T , ArrayXd &diffrho , double f )
{   ArrayXd T_T = (T > Tcrit).select(f*aK*(T-Tcrit) , 0);
    ArrayXd bracket = (P >= T_T).select(pow(P - T_T , perov_n-2) , -pow(T_T - P  , perov_n-2));
    diffrho = (T > Tcrit).select(perov_n * (perov_n -1)* perov_C * bracket * -f*aK, 0);

}

#endif 