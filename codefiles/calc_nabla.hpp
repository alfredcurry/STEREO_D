#ifndef _NABLA_CALC_H
#define _NABLA_CALC_H

#include <Eigen/Dense>
#include "StructStruct.hpp"
#include "EquationOfState.hpp"

/// @brief Class for calculing the log(T)/log(P) graidient
class nabla_calc 
{
    public:

        void nabla( struct P_P *planet); /// Calculates the gradient
        void diffnabla( struct P_P *planet0 , EoS *eos); /// Calculates the derivatives of the gradient w.r.t. P,r,T,L
        nabla_calc(int _I , int _J);
        void resize(int _J){/// Resizes
            J = _J;
        }
     
    private:

        int I , J;
        double Re_crit;
        void AB(struct P_P *planet , double &A , double &B , int j); /// Calculates coeffients for the super-adiabatic part of the gradient
        double nabnonadcalc(struct P_P *planet , const double &A , const double &B , const double &C , double Re_crit , int j ); /// Calculates the super-adiabatic part of the gradient
        double C_grav(struct P_P *planet , int j  );  /// Calculates coeffients associated with the gravitational part of the super-adiabatic gradient (Abe)
       
};
#endif