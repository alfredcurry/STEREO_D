#include "EquationOfState.hpp"
#include "EquationOfStateIdeal.hpp"
#include "2Dinterpolator.hpp"
#include <Eigen/Dense>

template <typename Real>
MeltFraction<Real>::MeltFraction(){}

template<class Real>
void EoS::compute_values( struct P_P *planet ){
    int *indexes;
    
    for(int j =0 ; j < planet->rho.size() ; j++){
        indexes = interp_rho_s.return_index(log10(planet->ym(j,0)),log10(planet->ym(j,2)));  
        planet->rho(j) = pow(10 , interp_rho_s.return_value(log10(planet->ym(j,0)),log10(planet->ym(j,2)),indexes[0],indexes[1]));
        planet->del_rho(j) = 1.0/planet->rho(j); //interp_del.return_value(log10(planet->ym(j,0)),planet->ym(j,2),indexes[0],indexes[1]);
        planet->Cp(j) = Ideal<double>::Cp(planet->ym(j,2),variables[1],variables[0]);
        planet->nu(j) = 1.0e-35; //For the ideal gas we assume a very low viscosity so that the mixing lenght regime is purely in the inviscid regime (which has no viscosity dependence)
        if( j == planet->rho.size()-1){
           // std::cout << planet->rho(j) << " rho " << Ideal<double>::EoSrho(planet->ym(j,0),planet->ym(j,2),1) << std::endl;
        }
    }
}

template<class Real>
void EoS::compute_derivs( struct P_P *planet ){
    int *indexes;
    double *grads;
    for(int j =0 ; j < planet->rho.size() ; j++){
        indexes = interp_rho_s.return_index(log10(planet->ym(j,0)),log10(planet->ym(j,2)));
        grads = interp_rho_s.return_grads(log10(planet->ym(j,0)),log10(planet->ym(j,2)),indexes[0],indexes[1]);
        planet->diffrhoP(j) = planet->rho(j)/planet->ym(j,0)*grads[0];
        planet->diffrhoT(j) = planet->rho(j)/planet->ym(j,2)*grads[1];
        
        //grads = interp_del.return_grads(log10(planet->ym(j,0)),planet->ym(j,2),indexes[0],indexes[1]);
        planet->diffdel_rhoP(j) = -1.0/(planet->rho(j)*planet->rho(j)) * planet->diffrhoP(j) ;//grads[0]/(P(j)*log(10));
        planet->diffdel_rhoT(j) = -1.0/(planet->rho(j)*planet->rho(j)) * planet->diffrhoT(j) ;//grads[1];
        planet->diffCpT(j) = 0;
        if( j == planet->rho.size()-1){
            std::cout << planet->diffrhoT(j) << " diffrho " << Ideal<double>::DiffrhoT(planet->ym(j,0),planet->ym(j,2),1) << std::endl;
        }
    }
}

EoS::EoS(double *_variables):
    visc({0})
{
    EoStype = "Ideal_Tab";
    variables.push_back(_variables[0]);
    variables.push_back(_variables[1]);
    std::string filebase , x_name , y_name , f_name;
    int grid_size[2];
    grid_size[0] = 1000;
    grid_size[1] = 500;

    x_name = "log_P";
    y_name = "log_T";
    filebase = "eos/Ideal";
    char mu_string[10];

    snprintf(mu_string,10 , "mu%.1f", variables[0]);
    filebase.append(mu_string);
    f_name = "log_rho";
    std::cout << "Preparing EoS...\n" << "From folder:\n" << filebase << std::endl;
    interp_rho_s.setup_from_coeff(filebase, x_name, y_name, f_name, grid_size);
    //f_name = "del";
    //interp_del.setup_from_coeff(filebase, x_name, y_name, f_name);
    std::cout << "EoS setup complete." << std::endl;
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