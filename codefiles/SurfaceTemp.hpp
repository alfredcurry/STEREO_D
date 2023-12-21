#ifndef _SURFACE_TEMP_H
#define _SURFACE_TEMP_H

#include <iostream>
#include <Eigen/Dense>
#include "PhysicalConsts.h"

template <class I>
class T_surface{
    public:
    
        static I BB_temp(const I L, const I R )
        {   
            return pow(L/(4.0*M_PI*pow(R,2)*Sig_Stefan),1.0/4);  
        }
        static I BB_diff_L(const I L, const I T)
        {
            return 0.25 * T / L;
        }
        static I BB_diff_R(const I R, const I T)
        {
            return -0.5 * T / R;
        }
        static I IRR_temp(const I Teff, const I Teq , const I tau , const I epsilon )
        {   
            if( Teff < 0){
               return Teq * ( 1.0 + 3.0/16 *  pow(Teff/Teq,4)/epsilon*(2.0/3 + tau));
            }else{ 
                return pow(3.0/4 * pow(Teff,4)/epsilon*(2.0/3 + tau) + pow(Teq,4) , 1.0/4);  
            }
        }
        static I IRR_diff_L(const I L, const I T , const I Teq , const I Teff , const I tau , const I epsilon)
        {   
            if( Teff < 0.1*Teq){
               return 12.0/16*pow(Teff/Teq,3)/epsilon*(2.0/3 + tau) * BB_diff_L( L, Teff);
            }else{ 
                return 3.0/4 * pow(Teff/T , 3) * (2.0 /3 + tau)/epsilon * BB_diff_L( L, Teff);
            }
        }
        static I IRR_diff_R(const I R, const I T , const I Teq , const I Teff , const I tau , const I epsilon)
        {
            if( Teff < 0.001*Teq){
               return 12.0/16*pow(Teff/Teq,3)/epsilon*(2.0/3 + tau) * BB_diff_R( R, Teff);
            }else{ 
                return 3.0/4 * pow(Teff/T , 3) * (2.0 /3 + tau)/epsilon * BB_diff_R( R , Teff);
            }
        }

        static I night_temp( const I &L , const I &R , const I &Tss){
            I Tn = pow( L / (2*M_PI * R * R * Sig_Stefan) , 0.25);
            if(Tn < Tss){
                return Tn;
            }else{
                return pow( L / (4*M_PI * R * R * Sig_Stefan) + 0.5*pow(Tss,4) , 0.25);
            }
        }
        static I diff_night_temp_dL( const I &Tn , const I &L , const I &R , const I &Tss){
            if(Tn < Tss){
                return 0.25*Tn/L;
            }else{
                return 0.25*pow(Tn,-3) / (4*M_PI * R * R * Sig_Stefan) ;
            }
        }

        static I diff_night_temp_dR( const I &Tn , const I &L , const I &R , const I &Tss){
            if(Tn < Tss){
                return -0.5*Tn/R;
            }else{
                return -0.5*pow(Tn,-3) * L / (4*M_PI * R* R * R * Sig_Stefan) ;
            }
        }
        static I z_coeffiencient( const I &T , const I &R , const I &m , const I &rho , const I k_cond , const I del_rho , const I &Cp , const I &visc , const I &DeltaT , const I &Pe ){
            return pow(visc*T/(pow(rho,3) * k_cond*k_cond *Cp * del_rho * G_Newt * pow(R,4) * m * DeltaT) * pow(4*M_PI * Pe , 2) , 1.0/3) / (4*M_PI);
        }
        static I param_bound_temp(const I &L , const I &z_coeff ,const I &Tss ,  const I &Tn){
            return L * z_coeff + 2.0/5 * Tss + 3.0/5*Tn;
        }
};
#endif
