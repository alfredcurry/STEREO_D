#include "EquationOfState.hpp"
#include "EquationOfStateIdeal.hpp"
#include <Eigen/Dense>
#include "Viscosity.hpp"

template <typename Real>
MeltFraction<Real>::MeltFraction(){}

template<class Real>
void EoS::compute_values( struct P_P *planet ){

    planet->rho = Ideal<Real>::EoSrho(planet->ym.col(0),planet->ym.col(2),variables[0]);
    planet->del_rho = Ideal<Real>::del(planet->ym.col(0),planet->ym.col(2),variables[0]) / planet->rho;
    planet->Cp = Ideal<Real>::Cp(planet->ym.col(2),variables[1],variables[0]);
    planet->nu = Eigen::ArrayXd::Constant(planet->rho.size(),1.0e-35);  //For the ideal gas we assume a very low viscosity so that the mixing lenght regime is purely in the inviscid regime (which has no viscosity dependence)
}

template<class Real>
void EoS::compute_derivs( struct P_P *planet ){
    planet->diffrhoP = Ideal<Real>::DiffrhoP(planet->ym.col(0),planet->ym.col(2),variables[0]);
    planet->diffrhoT = Ideal<Real>::DiffrhoT(planet->ym.col(0),planet->ym.col(2),variables[0]);
    planet->diffdel_rhoP = Ideal<Real>::DiffdelP(planet->ym.col(0),planet->ym.col(2),variables[0]) / planet->rho - planet->diffrhoP * planet->del_rho/planet->rho;
    planet->diffdel_rhoT = Ideal<Real>::DiffdelT(planet->ym.col(0),planet->ym.col(2),variables[0]) / planet->rho - planet->diffrhoT * planet->del_rho/planet->rho;
    planet->diffCpT = Ideal<Real>::CpDiffT(planet->ym.col(2));

}

EoS::EoS(double *_variables):
    visc({0})
{
    EoStype = "Ideal";
    variables.push_back(_variables[0]);
    variables.push_back(_variables[1]);

};

template <typename Real>
Viscosity<Real>::Viscosity(std::vector<Real> _variables)
{}

void TemporaryFunction ()
{
    Eigen::ArrayXd temp;
    double *a;
    EoS temp_EoS(a);
    P_P temp_P(0,0);
    temp_EoS.compute_values<Eigen::ArrayXd>(&temp_P);
    temp_EoS.compute_derivs<Eigen::ArrayXd>(&temp_P);    
}