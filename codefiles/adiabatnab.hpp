#ifndef _NAB_AD_H
#define _NAB_AD_H
#include "StructStruct.hpp"

template <class N> /// Class for calculating the adiabatic dlogT/dlogP gradient (nabla)
class Nab_Ad{

    public:
        static N nab(const N &P , const N &T , const N &del_rho , const N &Cp ){
            return P * del_rho/(T * Cp); 
        }
        static N DiffnabP(const N &P, const N &T, const N &del_rho , const N &Cp , const N & diffdel_rhoP , const N &diffCpP){
            //std::cout << P*diffdel_rhoP/(T*Cp) << " NABLA_DP " << del_rho/(T*Cp) << " TERMS " << - P*del_rho*diffCpP/(Cp*Cp*T) << std::endl;
            return P*diffdel_rhoP/(T*Cp) + del_rho/(T*Cp) - P*del_rho*diffCpP/(Cp*Cp*T);
        }
        static N DiffnabT(const N &P, const N &T , const N &del_rho , const N &Cp , const N &diffdel_rhoT , const N & diffCpT){
            
            return  P*diffdel_rhoT/(T*Cp) -P*del_rho/(T*T*Cp)  - P*del_rho*diffCpT/(Cp*Cp*T)  ;
        }
};

#endif //_NAB_AD_H