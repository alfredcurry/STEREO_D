#ifndef _Conductivity_H
#define _Conductivity_H

#include <iostream>
#include "PhysicalConsts.h"

template <class I> /// Class for calculating conductivity (assuming power-law behaviour)
class Conduct{
    public:
        static I k(const I rho, const I T , const double k_0, const double alpha , const double beta )
        {   
             return k_0 * pow(rho , alpha) * pow( T , beta) ;  
        }
        static I DiffkP(const I rho, const I diffrhoP, const I T , const double k_0, const double alpha , const double beta)
        { 
            return  k_0 * alpha* pow(rho, alpha - 1) * pow(T , beta) * diffrhoP ;
        }
        static I DiffkT(const I rho, const I del , const I T , const double k_0, const double alpha , const double beta)
        { 
            return k_0 * (beta* pow(rho, alpha ) * pow(T , beta-1) + alpha * pow(rho, alpha - 1 ) * -del * rho/T *pow(T, beta) );
        }
};
#endif