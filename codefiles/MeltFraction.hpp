#ifndef _MELT_FRAC_H
#define _MELT_FRAC_H

#include <Eigen/Dense>
#include <boost/math/interpolators/cubic_hermite.hpp>
/* 
in order to get the quintic b spline to work in the way I want I need to add the line 

    cardinal_cubic_hermite() = default;

to the "boost/math/interpolators/cubic_hermite.hpp" file after the constrcutor declarations in the cardinal_cubic_hermite class
*/
template <typename Real> ///Class for calculating the melt fraction, given T,P and some solius/liquidus curves
class MeltFraction
{
    private:
        Real Solidus(Real P);
        Real Liquidus(Real P);
   
        Eigen::Array<Real , Eigen::Dynamic , 1 >  NormTdP , NormTdT;
        std::vector<double> smooth_parameters; /// Non physical parameters for smooth approximations
        std::vector<double> analytical_parameters; /// Parameters for an analytical fit to solidus/liquidus

        boost::math::interpolators::cardinal_cubic_hermite<std::vector<Real>> Solidus_spline;
        boost::math::interpolators::cardinal_cubic_hermite<std::vector<Real>> Liquidus_spline;
        

    public:
        Real DiffTliqdP(Real P);
        Real DiffTsoldP(Real P);
        Eigen::Array<Real , Eigen::Dynamic , 1 > NormT , Tliq , Tsol , TsoldP , TliqdP , diff2phiPT , diff2phiTT  ; 
        Eigen::Array<Real , Eigen::Dynamic , 1 > Phi(Eigen::Array<Real , Eigen::Dynamic , 1 > P, Eigen::Array<Real , Eigen::Dynamic , 1 > T);
        Eigen::Array<Real , Eigen::Dynamic , 1 > DiffPhiT(Eigen::Array<Real , Eigen::Dynamic , 1 > T);
        Eigen::Array<Real , Eigen::Dynamic , 1 > DiffPhiP(Eigen::Array<Real , Eigen::Dynamic , 1 > P, Eigen::Array<Real , Eigen::Dynamic , 1 > T);
        Real DeltaS(Real DeltaV , Real phi , Real TsoldP , Real TliqdP); /// Calculates the specific entropy change of melting
        Real DiffDeltaSdP(Real DeltaS , Real DeltaV , Real DeltaVdT , Real diffphidT , Real TsoldP , Real TliqdP);
        Real DiffDeltaSdT(Real DeltaS , Real DeltaV , Real DeltaVdP , Real diffphidP , Real TsoldP , Real TliqdP);
        
        void SecondDerivs(Eigen::Array<Real , Eigen::Dynamic , 1 > P , Eigen::Array<Real , Eigen::Dynamic , 1 > T);
        MeltFraction();
};

#endif