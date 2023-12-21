#ifndef _B_C_H
#define _B_C_H

#include <vector>
#include <Eigen/Dense>
#include "2Dinterpolator.hpp"
#include "fitted_BC.hpp"
#include <string>
#include "StructStruct.hpp"

/// @brief Class for calculating the temperature boundary condition
class B_C{
    public:
        template <class Real> /// Calculates the temperature at the edge of the grid
        Real T_surf(struct P_P *planet , Real &L);
        template <class Real> /// Calcualates derivatives (w.r.t. P,r,T,L) of the temperature at the edge of the grid
        void T_surf_derivs(struct P_P *planet );
        void read_table(double &P0 , double &Tss); /// If required, reads in data from T_surf parameterisation
        B_C();
        char BCtype;
        
    private:
        Interpolators::BicubicInterpolator<double> interp_T0;
        Fitted_Layer fit;
};


#endif