#ifndef _EoS_CALC_H
#define _EoS_CALC_H

#include <vector>
#include <Eigen/Dense>
//#include "EquationOfStateIdeal.hpp"
#include "2Dinterpolator.hpp"
#include <string>
#include "StructStruct.hpp"
#include "MeltFraction.hpp"
#include "Viscosity.hpp"
#include <boost/math/interpolators/pchip.hpp>

/// Class for evaluation the equation of state. EoS type determined at compile time, with certain parameter tweaks available at run-time
class EoS{
    public:
        template <class Real> /// Computes the values of the EoS at each point in the planet
        void compute_values(struct P_P *planet );
        template <class Real>  /// Computes the derivatives the EoS values at each point in the planet
        void compute_derivs(struct P_P *planet );
        template <class Real>  /// Computes the values of the density only
        void compute_density_only( struct P_P *planet );
        EoS(double *_variables);
        std::string EoStype;
        MeltFraction<double> melt;
        Viscosity<double> visc;

    private:
        int J;
        std::vector<double> variables;
        Interpolators::BicubicInterpolator<double> interp_rho_m , interp_del_rho_m , interp_Cp_m , interp_rho_s , interp_del_rho_s , interp_Cp_s , interp_rho_c , interp_del_rho_c , interp_Cp_c , interp_rho_cl , interp_del_rho_cl , interp_Cp_cl;
        boost::math::interpolators::pchip<std::vector<double>> Core_melt_spline;
};


#endif