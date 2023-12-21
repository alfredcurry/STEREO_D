#ifndef _EQUATION_OF_STATE_H
#define _EQUATION_OF_STATE_H

#include <iostream>
#include "PhysicalConsts.h"

template <class I>
class Ideal{
    public:
        static I EoSrho(const I P, const I T , const double mu_m)
        {   
             return P / (T * boltz_R/mu_m) ;  
        }
        static I DiffrhoP(const I P, const I T , const double mu_m)
        { 
            return  1 / (T*boltz_R/mu_m) ;
        }
        static I DiffrhoT(const I P, const I T , const double mu_m)
        { 
            return - P / (pow(T,2)*boltz_R/mu_m) ;
        }
        static I del(const I P, const I T , const double mu_m)
        {
            return pow(P,0);
        }
         static I DiffdelP(const I P, const I T , const double mu_m)
        {
            return 0*P;
        }
        static I DiffdelT(const I P, const I T , const double mu_m)
        {
            return 0*P;
        }
        static I Diff2delPT(const I P, const I T , const double mu_m)
        {
            return 0*P;
        }
        static I Diff2rhoTT(const I P, const I T ,  const double mu_m)
        { 
            return 2*P / (pow(T,3)*boltz_R/mu_m) ;
        }
        static I Diff2rhoPT(const I P, const I T, const double mu_m)
        { 
            return -1/(pow(T,2)*boltz_R/mu_m) ;
        }
        static I Cp(const I T , const double n , const double mu_m )
        {
            double gamma = 1+1/n ; 
            return gamma/(gamma-1)*boltz_R/mu_m * pow(T,0);
        }
        static I CpDiffT(const I T )
        {
            return 0 * T ;
        }
};
#endif 