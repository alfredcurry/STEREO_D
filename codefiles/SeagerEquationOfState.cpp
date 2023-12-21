#include "EquationOfState.hpp"
#include "EquationOfStateSea.hpp"
#include <Eigen/Dense>
#include "Viscosity.hpp"

template <typename Real>
MeltFraction<Real>::MeltFraction(){}

template<class Real>
void EoS::compute_values( struct P_P *planet ){
    EoSrhosea(planet->ym.col(0),planet->ym.col(2), planet->rho ,variables[0]);
    DiffrhoTsea(planet->ym.col(0),planet->ym.col(2), planet->del_rho ,variables[0]);
    planet->del_rho = planet->del_rho*planet->ym.col(2)/(planet->rho*planet->rho);
    planet->Cp = 1200;
}

template<class Real>
void EoS::compute_derivs( struct P_P *planet ){
    DiffrhoPsea(planet->ym.col(0),planet->ym.col(1),planet->diffrhoP,variables[0]);
    DiffrhoTsea(planet->ym.col(0),planet->ym.col(1),planet->diffrhoT,variables[0]);
    Diff2rhoPTsea(planet->ym.col(0),planet->ym.col(1),planet->diffdel_rhoP,variables[0]);
    planet->diffdel_rhoP = -planet->ym.col(1)/(planet->rho*planet->rho) * planet->diffrhoP*planet->diffrhoT + planet->ym.col(1)/planet->rho * planet->diffdel_rhoP;
    Diff2rhoTTsea(planet->ym.col(0),planet->ym.col(1),planet->diffdel_rhoT,variables[0]);
    planet->diffdel_rhoT = planet->diffrhoT/planet->rho - planet->ym.col(1)/(planet->rho*planet->rho) *planet->diffrhoT*planet->diffrhoT + planet->ym.col(1)/planet->rho *planet->diffdel_rhoT;
    
    planet->diffCpT = 0*planet->ym.col(1);

}

EoS::EoS(double *_variables):
    visc({0})
{
    EoStype = "Seager";
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