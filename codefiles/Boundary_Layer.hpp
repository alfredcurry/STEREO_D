#ifndef _BL_H
#define _BL_H

#include <vector>
#include "PhysicalConsts.h"

template <typename Real> /// Class used for parameterisation of an unresolved boundary layer at the edge of the domain
class Boundary_Layer
{
    private: 
        Real Pe_crit = 3.0;

    public:
        Real DeltaT(Real &L , Real &R , const Real &visc , Real &T , Real &del_rho , Real &rho , Real &m , Real &k , Real &Cp){
            return pow(L , 0.75)*pow( Pe_crit*Pe_crit*visc*T/(4*M_PI * pow(R,4) * del_rho*pow(rho,3)*G_Newt*m*k*k*Cp) , 0.25);
        }
        Real DeltaP(Real &rho , Real &m , Real &DeltaT, Real &L , Real &k){
            return 4*M_PI*rho*G_Newt*m*DeltaT*k/L;
        }
        Real Deltaz(Real &DeltaT , Real &k , Real &L , Real  &R){
            return 4*M_PI*R*R * k * DeltaT/L ;
        }

        static Real DiffDeltaT_P( Real &DeltaT , Real &rho , Real &drhoP , Real &del_rho  , Real &ddel_rhoP , Real &Cp ,  Real &dCpP , const Real &visc , const Real &dviscP ){
            return DeltaT*( -0.75/rho * drhoP - 0.25/del_rho * ddel_rhoP - 0.25/Cp * dCpP + 0.25/visc * dviscP ); 
        }
        static Real DiffDeltaT_R( Real &DeltaT , Real &R){
            return - DeltaT/R;
        }
        static Real DiffDeltaT_T( Real &DeltaT , Real &rho , Real &drhoT , Real &del_rho  , Real &ddel_rhoT , Real &Cp , Real &dCpT , const Real &visc , const Real &dviscT , Real &T){
            return DeltaT*( -0.75/rho * drhoT - 0.25/del_rho * ddel_rhoT - 0.25/Cp * dCpT + 0.25/visc * dviscT + 0.25/T); 
        }
        static Real DiffDeltaT_L(Real &DeltaT , Real &L){
            return 0.75*DeltaT/L;
        }
        
        static Real DiffDeltaP_P(Real &DeltaP , Real &rho , Real &drhoP , Real &DeltaT , Real &dDeltaT_P){
            return DeltaP * ( drhoP/rho + dDeltaT_P/DeltaT );
        }
        static Real DiffDeltaP_R(Real &DeltaP , Real &DeltaT , Real &dDeltaT_R){
            return DeltaP * dDeltaT_R/DeltaT ;
        }
        static Real DiffDeltaP_T(Real &DeltaP , Real &rho , Real &drhoT , Real &DeltaT , Real &dDeltaT_T){
            return DeltaP * ( drhoT/rho + dDeltaT_T/DeltaT );
        }
        static Real DiffDeltaP_L(Real &DeltaP , Real &L , Real &DeltaT , Real &dDeltaT_L){
            return DeltaP * (-1.0/L + dDeltaT_L/DeltaT );
        }

        static Real DiffDeltaz_P( Real &Deltaz , Real &DeltaT , Real &dDeltaT_P){
            return Deltaz/DeltaT * dDeltaT_P;
        }
        static Real DiffDeltaz_R( Real &Deltaz , Real &R , Real &DeltaT , Real &dDeltaT_R ){
            return 2*Deltaz/R + Deltaz/DeltaT * dDeltaT_R;
        }
        static Real DiffDeltaz_T( Real &Deltaz , Real &DeltaT , Real &dDeltaT_T){
            return Deltaz/DeltaT * dDeltaT_T;
        }
        static Real DiffDeltaz_L( Real &Deltaz , Real &L , Real &DeltaT , Real &dDeltaT_L){
            return - Deltaz/L + Deltaz/DeltaT * dDeltaT_L;
        }
};

#endif