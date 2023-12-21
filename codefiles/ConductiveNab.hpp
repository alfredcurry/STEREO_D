#ifndef _NAB_cond_H
#define _NAB_cond_H

#include <iostream>
#include "StructStruct.hpp"

template <class N> /// Class for calculating the logTlogP temperature gradient (nabla) for a conductive (or radiative) structure
class CondNab
{
    public:
        static N nab(const N &L , const N &K){
            return L / K;
        }
        N diffnabP(const N &L , const N &K , const N &P , const N &rho , const N &diffrhoP , const N k , const N &diffkP){
            return L/K * (1.0/P - 1.0/rho*diffrhoP - 1/k*diffkP);
        }
        N diffnabT(const N &L , const N &T , const N &K , const N &rho , const N &diffrhoT , const N k , const N &diffkT ){
            return L/K * (-1.0/T + 1.0/rho*diffrhoT -1/k*diffkT);
        }
        N diffnabL(const N &K ){
            return 1.0/K; //although bear with if there's L dependece on opacity/conductivity
        }
};

#endif