#include "EquationOfState.hpp"
#include "2Dinterpolator.hpp"
#include <Eigen/Dense>

template <typename Real>
MeltFraction<Real>::MeltFraction(){}

template<class Real>
void EoS::compute_values( struct P_P *planet ){
    int *indexes;

    for(int j =0 ; j < planet->rho.size() ; j++){
        indexes = interp_rho_s.return_index(log10((planet->ym(j,0)/1e9)),log10(planet->ym(j,2)));  

        planet->rho(j) = pow(10 , interp_rho_s.return_value(log10((planet->ym(j,0)/1e9)),log10(planet->ym(j,2)),indexes[0],indexes[1]));

        planet->del_rho(j) = pow(10 , interp_del_rho_s.return_value(log10((planet->ym(j,0)/1e9)),log10(planet->ym(j,2)),indexes[0],indexes[1]));
        planet->Cp(j) = pow(10 , interp_Cp_s.return_value(log10((planet->ym(j,0)/1e9)),log10(planet->ym(j,2)),indexes[0],indexes[1]));
        if( j == planet->rho.size()-1){
           // std::cout << rho(j) << " rho " << Ideal<double>::EoSrho((planet->ym(j,0)/1e9),planet->ym(j,2),1) << std::endl;
        }
    }
}

template<class Real>
void EoS::compute_derivs( struct P_P *planet ){
    int *indexes;
    double *grads;
    for(int j =0 ; j < planet->rho.size() ; j++){

        indexes = interp_rho_s.return_index(log10((planet->ym(j,0)/1e9)),log10(planet->ym(j,2)));

        grads = interp_rho_s.return_grads(log10((planet->ym(j,0)/1e9)),log10(planet->ym(j,2)),indexes[0],indexes[1]);
        planet->diffrhoP(j) = planet->rho(j)/planet->ym(j,0)*grads[0];
        planet->diffrhoT(j) = planet->rho(j)/planet->ym(j,2)*grads[1];
        
        grads = interp_del_rho_s.return_grads(log10((planet->ym(j,0)/1e9)),log10(planet->ym(j,2)),indexes[0],indexes[1]);
        planet->diffdel_rhoP(j) = planet->del_rho(j)/planet->ym(j,0)* grads[0] ;
        planet->diffdel_rhoT(j) = planet->del_rho(j)/planet->ym(j,2)* grads[1] ;
        
        grads = interp_Cp_s.return_grads(log10((planet->ym(j,0)/1e9)),log10(planet->ym(j,2)),indexes[0],indexes[1]);
        planet->diffCpP(j) = planet->Cp(j)/planet->ym(j,0)* grads[0] ;
        planet->diffCpT(j) = planet->Cp(j)/planet->ym(j,2)* grads[1] ;
      
    }
}

EoS::EoS(double *_variables):
    visc({0})
{
    EoStype = "Stixrude";
    
    std::string filebase , x_name , y_name , f_name;
    int grid_size[2];
    grid_size[0] = 250;
    grid_size[1] = 250;

    x_name = "logP";
    y_name = "logT";
    filebase = "eos/Stixrude";
   
    std::cout << "Preparing EoS...\n" << "From folder:\n" << filebase << std::endl;

    f_name = "logrho";
    std::cout << "Preparing density..." << std::endl;
    interp_rho_s.setup_from_coeff(filebase, x_name, y_name, f_name, grid_size);
    
    f_name = "logdel_rho";
    std::cout << "Preparing delta..." << std::endl;
    interp_del_rho_s.setup_from_coeff(filebase, x_name, y_name, f_name, grid_size);
    
    f_name = "logCp";
    std::cout << "Preparing heat capacity..." << std::endl;
    interp_Cp_s.setup_from_coeff(filebase, x_name, y_name, f_name, grid_size);

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