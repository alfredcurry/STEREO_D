#ifndef _IRON_CORE_H
#define _IRON_CORE_H

template <class Real> /// Class for calculating the heat input from an UNRESOLVED iron core 
class Iron_core {

    private:
        Real m , c = 880, T_adjust = 1.47;
    public:
        Real Lcore( Real T , Real T0 , Real dt){
            return - m * c * T_adjust * (T- T0) / dt ;
        }
        Real dLcore_dT( Real dt ){
            return - m * c * T_adjust / dt ;
        }
        Iron_core(Real _m):
        m(_m) {}

};

#endif