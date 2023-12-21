#include "Viscosity.hpp"
#include <vector>
#include "PhysicalConsts.h"

/**
 * @brief Implementation of the Viscosity class for viscosity for molten/solid rock. Combines variability with melt fraction with an Arrhenius law for the solid viscosity, and constant viscosity for the melt. The code does NOT need the dervatives, so these may be out-of-date.
 * 
 */


template <typename Real>
Real Tanh_step(Real &x , Real &x0 , Real &k){
    return 0.5 + 0.5*tanh(k*(x-x0));
}
template <typename Real>
Real Tanh_step_deriv(Real &x , Real &x0 , Real &k){
    return 0.5*k*(1.0-tanh(k*(x-x0))*tanh(k*(x-x0)));
}

template <typename Real> /// Arrhenius law type viscosity function
Real Arrhenius_Viscosity( Real &P ,Real &T , Real &visc_0 , Real&T0 , Real &V0 , Real&E0 ){

    return visc_0 * exp( ((E0 + P*V0)/T - E0/T0) / R_gas);
}
template <typename Real>
Real Diff_Arrhenius_Viscosity_dT( const Real &visc ,Real &P , Real &T , Real &V0 , Real&E0 ){

    return -visc * (E0 + P*V0)/(R_gas*T*T);
}
template <typename Real>
Real Diff_Arrhenius_Viscosity_dP( const Real &visc , Real &T , Real &V0 ){

    return visc * V0/(R_gas*T);
}

// variables are phi_c , alpha , visc_l , visc_s(0) [ Pa s] , (T0 , V0 [m^3mol-1], E0 [ Jmol-1])*/
template <typename Real>
Real Viscosity<Real>::Visc(Real &phi , Real &P , Real &T){
    //if no melting T refers to the actual T. If there is melting T refers to T at phi=0 for that pressure.
    Real visc_s;
    if(variables[9] < 0.0){
        visc_s = variables[3];
    }else{
        visc_s = Arrhenius_Viscosity(P,T, variables[3] , variables[4] , variables[5] , variables[6]) ;
    }
    
    if( phi == 0){
        return visc_s;
    }else
    if(phi == 1){
        return variables[2];
    }else/*{
        return 1.0/( exp(variables[1]*phi)/Arrhenius_Viscosity(P,T, variables[3] , variables[4] , variables[5] , variables[6]) + Tanh_step(phi,variables[0], variables[7]) * ( pow((phi-variables[0])*(phi-variables[0]),5.0/4)/pow(1-variables[0],2.5)/variables[2] + 1.0/variables[8]));
    }*/
    if(phi < variables[0]){
        return visc_s*exp(-variables[1]*phi);
    }else
    {
        return 1.0/( exp(variables[1]*phi)/visc_s + pow((phi-variables[0])/(1-variables[0]),2.5)/variables[2]);
    }
}

/* THE DIFFS ARE UNUSED SO OUT OF DATE !*/

template <typename Real>
Real Viscosity<Real>::DiffVisc_dT(Real &phi , const Real &visc , Real &P , Real &T , Real &dphiT ){ /* dont need */
    if( phi == 0){
        return Diff_Arrhenius_Viscosity_dT( visc , P , T , variables[5] , variables[6] );
    }else
    if(phi == 1){
        return 0;
    }else/*{
        Real visc_s = Arrhenius_Viscosity(P,T, variables[3] , variables[4] , variables[5] , variables[6]);
        return -(variables[1]*exp(variables[1]*phi)/visc_s*dphiT - exp(variables[1]*phi)/(visc_s*visc_s)*Diff_Arrhenius_Viscosity_dT( visc , P , T , variables[5] , variables[6]) 
            + dphiT * (Tanh_step_deriv(phi,variables[0], variables[7])* (pow((phi-variables[0])*(phi-variables[0]),5.0/4)/pow(1-variables[0],2.5)/variables[2] + 1.0/variables[8]) + Tanh_step(phi,variables[0], variables[7]) * 2.5 * (phi-variables[0]) * pow((phi-variables[0])*(phi-variables[0]),0.25)/pow(1-variables[0],2.5)/variables[2])
                )/pow( exp(variables[1]*phi)/visc_s + Tanh_step(phi,variables[0], variables[7]) * (pow((phi-variables[0])*(phi-variables[0]),5.0/4)/pow(1-variables[0],2.5)/variables[2] + 1.0/variables[8]), 2);

    }*/
    if(phi < variables[0]){
        return -variables[1]*Arrhenius_Viscosity(P,T, variables[3] , variables[4] , variables[5] , variables[6])*exp(-variables[1]*phi)*dphiT ;
    }else
    {   
        double visc_s = Arrhenius_Viscosity(P,T, variables[3] , variables[4] , variables[5] , variables[6]);
        return -(variables[1]*exp(variables[1]*phi)/visc_s*dphiT + 2.5*pow(phi-variables[0], 1.5)/pow(1-variables[0],2.5)/variables[2]*dphiT )/pow( exp(variables[1]*phi)/visc_s + pow((phi-variables[0])/(1-variables[0]),2.5)/variables[2] , 2) ;
    }
}
template <typename Real>
Real Viscosity<Real>::DiffVisc_dP(Real &phi , const Real &visc , Real &P , Real &T , Real &dphiP , Real &dT_phi0dP ){ /* dont need unless have BC */
    if( phi == 0){
        return Diff_Arrhenius_Viscosity_dP( visc , T , variables[5] );
    }else
    if(phi == 1){
        return 0;
    }else/*{
        Real visc_s = Arrhenius_Viscosity(P,T, variables[3] , variables[4] , variables[5] , variables[6]);
        return -(variables[1]*exp(variables[1]*phi)/visc_s*dphiP - exp(variables[1]*phi)/(visc_s*visc_s)*Diff_Arrhenius_Viscosity_dP( visc ,  T , variables[5] ) 
            + dphiP * (Tanh_step_deriv(phi,variables[0], variables[7]) * (pow((phi-variables[0])*(phi-variables[0]),5.0/4)/pow(1-variables[0],2.5)/variables[2] + 1.0/variables[8]) + Tanh_step(phi,variables[0], variables[7]) * 2.5 * (phi-variables[0]) * pow((phi-variables[0])*(phi-variables[0]),0.25)/pow(1-variables[0],2.5)/variables[2])
                )/pow( exp(variables[1]*phi)/visc_s + Tanh_step(phi,variables[0], variables[7]) * (pow((phi-variables[0])*(phi-variables[0]),5.0/4)/pow(1-variables[0],2.5)/variables[2] + 1.0/variables[8]), 2);

    }*/
    if(phi < variables[0]){
        return -variables[1]*Arrhenius_Viscosity(P,T, variables[3] , variables[4] , variables[5] , variables[6])*exp(-variables[1]*phi)*dphiP + Diff_Arrhenius_Viscosity_dP( visc , T , variables[5] )*exp(-variables[1]*phi) + Diff_Arrhenius_Viscosity_dT( visc , P , T , variables[5] , variables[6])*exp(-variables[1]*phi)*dT_phi0dP;
    }else
    {   
        double visc_s = Arrhenius_Viscosity(P,T, variables[3] , variables[4] , variables[5] , variables[6]);
        return -(variables[1]*exp(variables[1]*phi)/visc_s*dphiP - exp(variables[1]*phi)/(visc_s*visc_s)*Diff_Arrhenius_Viscosity_dP( visc , T , variables[5]) - exp(variables[1]*phi)/(visc_s*visc_s)*Diff_Arrhenius_Viscosity_dT( visc , P , T , variables[5] , variables[6])*dT_phi0dP + 2.5*pow(phi-variables[0], 1.5)/pow(1-variables[0],2.5)/variables[2] *dphiP )/pow( exp(variables[1]*phi)/visc_s + pow((phi-variables[0])/(1-variables[0]),2.5)/variables[2] , 2) ;
    }
}

template <typename Real>
Viscosity<Real>::Viscosity(std::vector<Real> _variables): /// Constructor just reads in the variables
    variables(_variables)
{
        std::cout << "Viscosity set up." << std::endl;
}