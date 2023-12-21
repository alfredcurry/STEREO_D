#include "EquationOfState.hpp"
#include "MeltFraction.hpp"
#include "MeltFunctions.hpp"
#include <Eigen/Dense>
#include "StructStruct.hpp"
#include "Viscosity.hpp"
//#include "MeltViscosity.hpp"
#include "MeltViscosity_Arrhenius.hpp"
#include "Vector_reader.h"
#include <boost/math/interpolators/pchip.hpp>

template<class Real>
void EoS::compute_values( struct P_P *planet ){
    int *indexes;
    J = planet->rho.size() ;
    planet->phi.tail(J-planet->Jcore) = melt.Phi(planet->ym.block(planet->Jcore,0,J-planet->Jcore,1) , planet->ym.block(planet->Jcore,2,J-planet->Jcore,1));
    
    if(planet->Jcore > 0){
        planet->phi.head(planet->Jcore) = 0;
    }
    planet->diffphiT.tail(J-planet->Jcore) = melt.DiffPhiT(planet->ym.block(planet->Jcore,2,J-planet->Jcore,1));
    //std::cout << planet->ym.col(0) << " TP " << planet->ym.col(2) <<std::endl;
    planet->diffphiP.tail(J-planet->Jcore) = melt.DiffPhiP(planet->ym.block(planet->Jcore,0,J-planet->Jcore,1) , planet->ym.block(planet->Jcore,2,J-planet->Jcore,1)); 
    //planet->phi = Eigen::ArrayXd::Constant(planet->nabla.size(),1.0);
    
    double one = 1;
    planet->visc_m = visc.Visc(one,one,one);

    if( planet->Jcore > 0){
        for(int j = 0 ; j < planet->Jcore ; j++){
            indexes = interp_rho_c.return_index(planet->ym(j,0)/1e9,planet->ym(j,2));  

            planet->rho(j) = interp_rho_c.return_value(planet->ym(j,0)/1e9,planet->ym(j,2),indexes[0],indexes[1]);
            planet->del_rho(j) = interp_del_rho_c.return_value(planet->ym(j,0)/1e9,planet->ym(j,2),indexes[0],indexes[1]);
            planet->Cp(j) = 600;            
            planet->phi(j) = 0.0;
        }
    }
    for(int j = planet->Jcore ; j < J ; j++){
        if(planet->phi(j) < 1){ //If any solid need to work out solid values
            
            indexes = interp_rho_s.return_index(planet->ym(j,0)/1e9,planet->ym(j,2));  

            planet->rho_s(j) = interp_rho_s.return_value(planet->ym(j,0)/1e9,planet->ym(j,2),indexes[0],indexes[1]);
            planet->del_rho_s(j) = interp_del_rho_s.return_value(planet->ym(j,0)/1e9,planet->ym(j,2),indexes[0],indexes[1]);
            planet->Cp_s(j) = 1500;
            
        }
        if(planet->phi(j)>0){ //If any melt need to work out melt values
    
            indexes = interp_rho_m.return_index(planet->ym(j,0)/1e9,planet->ym(j,2));  
    
            planet->rho_m(j) = interp_rho_m.return_value(planet->ym(j,0)/1e9,planet->ym(j,2),indexes[0],indexes[1]);
            planet->del_rho_m(j) = interp_del_rho_m.return_value(planet->ym(j,0)/1e9,planet->ym(j,2),indexes[0],indexes[1]);
            planet->Cp_m(j) = 1500;
        }
        if(planet->phi(j)==0){
            planet->rho(j) = planet->rho_s(j);
            planet->del_rho(j) = planet->del_rho_s(j);
            planet->Cp(j) = planet->Cp_s(j);
            
            planet->del_rho_avg(j) = planet->del_rho(j);
            planet->Cp_avg(j) = planet->Cp(j);

            planet->nu(j) = std::min(visc.Visc(planet->phi(j),planet->ym(j,0),planet->ym(j,2))/planet->rho(j) , 1e35);

        }else 
        if(planet->phi(j)==1){
            planet->rho(j) = planet->rho_m(j);
            planet->del_rho(j) = planet->del_rho_m(j);
            planet->Cp(j) = planet->Cp_m(j);

            planet->del_rho_avg(j) = planet->del_rho(j);
            planet->Cp_avg(j) = planet->Cp(j);
            planet->nu(j) = visc.Visc(one,one,one)/planet->rho(j);
            
        }else{    

            planet->DeltaV(j) = 1.0/planet->rho_m(j) - 1.0/planet->rho_s(j) ;

            planet->DeltaS(j) = 4e5/planet->ym(j,2);

            if( planet->DeltaS(j) < 0){
                std::cout << "Delta S odd" << std::endl;
                std::cout << " P " << planet->ym(j,0)/1e9 << " T " << planet->ym(j,2) << " melt " << planet->rho_m(j) << " solid " << planet->rho_s(j) << " DV " << planet->DeltaV(j) << std::endl;
            }
            
            planet->rho(j) = 1.0 / (planet->phi(j)/planet->rho_m(j) + (1-planet->phi(j))/planet->rho_s(j));
            planet->del_rho_avg(j) = planet->phi(j)*planet->del_rho_m(j) + (1-planet->phi(j))*planet->del_rho_s(j);
            planet->del_rho(j) = planet->del_rho_avg(j) + planet->DeltaV(j)*planet->ym(j,2)*planet->diffphiT(j) ;     
            planet->Cp_avg(j) = planet->phi(j)*planet->Cp_m(j) + (1-planet->phi(j))*planet->Cp_s(j);
            planet->Cp(j) = planet->Cp_avg(j) + planet->ym(j,2)*planet->DeltaS(j)*planet->diffphiT(j) ;

            planet->nu(j) = std::min(visc.Visc(planet->phi(j),planet->ym(j,0),melt.Tsol(j-planet->Jcore))/planet->rho(j) , 1e35);
            
        }

    }
  
}

template<class Real>
void EoS::compute_derivs( struct P_P *planet ){
    int *indexes;
    double *grads;
    double diffrhoP_s , diffrhoT_s , diffrhoP_m , diffrhoT_m , diffdel_rhoP_s , diffdel_rhoT_s , diffdel_rhoP_m , diffdel_rhoT_m , diffCpP_s , diffCpT_s , diffCpP_m , diffCpT_m , diffDeltaVdP , diffDeltaVdT; 
    
    melt.SecondDerivs(planet->ym.block(planet->Jcore,0,J-planet->Jcore,1) , planet->ym.block(planet->Jcore,2,J-planet->Jcore,1));

    if( planet->Jcore > 0){
        for(int j =0 ; j < planet->Jcore ; j++){
            
            indexes = interp_rho_c.return_index(planet->ym(j,0)/1e9,planet->ym(j,2));

            grads = interp_rho_c.return_grads(planet->ym(j,0)/1e9,planet->ym(j,2),indexes[0],indexes[1]);      
            planet->diffrhoP(j) = grads[0]/1e9;
            planet->diffrhoT(j) = grads[1];
        
            grads = interp_del_rho_c.return_grads(planet->ym(j,0)/1e9,planet->ym(j,2),indexes[0],indexes[1]);
            planet->diffdel_rhoP(j) = grads[0]/1e9 ;
            planet->diffdel_rhoT(j) = grads[1] ;

            planet->diffCpP(j) = 0.0;
            planet->diffCpT(j) = 0.0;
        }
    }
    for(int j = planet->Jcore ; j < J; j++){
        if(planet->phi(j) < 1){ //If any solid need to work out solid values
            indexes = interp_rho_s.return_index(planet->ym(j,0)/1e9,planet->ym(j,2));

            grads = interp_rho_s.return_grads(planet->ym(j,0)/1e9,planet->ym(j,2),indexes[0],indexes[1]);
            
            diffrhoP_s = grads[0]/1e9;
            diffrhoT_s = grads[1];

            grads = interp_del_rho_s.return_grads(planet->ym(j,0)/1e9,planet->ym(j,2),indexes[0],indexes[1]);
            diffdel_rhoP_s = grads[0]/1e9 ;
            diffdel_rhoT_s = grads[1] ;

            diffCpP_s = 0.0;
            diffCpT_s = 0.0;

        }
        if(planet->phi(j)>0){ //If any melt need to work out melt values
            
            indexes = interp_rho_m.return_index(planet->ym(j,0)/1e9,planet->ym(j,2));  

            grads = interp_rho_m.return_grads(planet->ym(j,0)/1e9,planet->ym(j,2),indexes[0],indexes[1]);
            diffrhoP_m = grads[0]/1e9;
            diffrhoT_m = grads[1];
        
            grads = interp_del_rho_m.return_grads(planet->ym(j,0)/1e9,planet->ym(j,2),indexes[0],indexes[1]);
            diffdel_rhoP_m = grads[0]/1e9 ;
            diffdel_rhoT_m = grads[1] ;
            
            diffCpP_m = 0.0;
            diffCpT_m = 0.0;
        
        }
        if(planet->phi(j)==0){

            planet->diffrhoP(j) = diffrhoP_s;
            planet->diffrhoT(j) = diffrhoT_s;
            planet->diffdel_rhoP(j) = diffdel_rhoP_s;
            planet->diffdel_rhoT(j) = diffdel_rhoT_s;
            planet->diffCpP(j) = diffCpP_s;
            planet->diffCpT(j) = diffCpT_s;
            
        }else 
        if(planet->phi(j)==1){

            planet->diffrhoP(j) = diffrhoP_m;
            planet->diffrhoT(j) = diffrhoT_m;
            planet->diffdel_rhoP(j) = diffdel_rhoP_m;
            planet->diffdel_rhoT(j) = diffdel_rhoT_m;
            planet->diffCpP(j) = diffCpP_m;
            planet->diffCpT(j) = diffCpT_m;

        }else{

            diffDeltaVdP = diffrhoP_s/(planet->rho_s(j)*planet->rho_s(j)) - diffrhoP_m/(planet->rho_m(j)*planet->rho_m(j));
            diffDeltaVdT = diffrhoT_s/(planet->rho_s(j)*planet->rho_s(j)) - diffrhoT_m/(planet->rho_m(j)*planet->rho_m(j));
            
            planet->diffrhoP(j) = -planet->rho(j)*planet->rho(j) * ( planet->diffphiP(j)*planet->DeltaV(j)
                - planet->phi(j)/(planet->rho_m(j)*planet->rho_m(j)) * diffrhoP_m 
                    - (1-planet->phi(j))/(planet->rho_s(j)*planet->rho_s(j))*diffrhoP_s );
            planet->diffrhoT(j) = -planet->rho(j)*planet->rho(j) * ( planet->diffphiT(j)*planet->DeltaV(j)
                - planet->phi(j)/(planet->rho_m(j)*planet->rho_m(j)) * diffrhoT_m 
                    - (1-planet->phi(j))/(planet->rho_s(j)*planet->rho_s(j))*diffrhoT_s );
            
            planet->diffdel_rhoP(j) = planet->diffphiP(j)*(planet->del_rho_m(j)-planet->del_rho_s(j)) + planet->phi(j)*diffdel_rhoP_m + (1-planet->phi(j))*diffdel_rhoP_s 
                + planet->DeltaV(j) * planet->ym(j,2) * melt.diff2phiPT(j-planet->Jcore) 
                    + diffDeltaVdP * planet->ym(j,2)*planet->diffphiT(j);
            
            planet->diffdel_rhoT(j) = planet->diffphiT(j)*(planet->del_rho_m(j)-planet->del_rho_s(j)) + planet->phi(j)*diffdel_rhoT_m + (1-planet->phi(j))*diffdel_rhoT_s 
                + planet->DeltaV(j) * planet->diffphiT(j) + planet->DeltaV(j) * planet->ym(j,2) * melt.diff2phiTT(j-planet->Jcore) 
                    + diffDeltaVdT * planet->ym(j,2)*planet->diffphiT(j);
            
            planet->diffCpP(j) = planet->phi(j)*diffCpP_m + (1-planet->phi(j))*diffCpP_s + (planet->Cp_m(j)-planet->Cp_s(j))*planet->diffphiP(j) + planet->DeltaS(j)*planet->ym(j,2)*melt.diff2phiPT(j-planet->Jcore) ; 
            planet->diffCpT(j) = planet->phi(j)*diffCpT_m + (1-planet->phi(j))*diffCpT_s + (planet->Cp_m(j)-planet->Cp_s(j))*planet->diffphiT(j) + planet->DeltaS(j)*planet->ym(j,2)*melt.diff2phiTT(j-planet->Jcore) ;
           
        }
    }

    planet->diffrhoP = planet->diffrhoP ;
    planet->diffCpP = planet->diffCpP ;
    planet->diffdel_rhoP = planet->diffdel_rhoP ;

    //planet->dviscdphi = visc.DiffVisc(planet->phi(J-1));
//    planet->dviscdP = visc.DiffVisc_dP(planet->phi(J-1) , planet->nu(J-1)*planet->rho(J-1) , planet->ym(J-1,0) , planet->ym(J-1,2) , planet->diffphiP(J-1));
  //  planet->dviscdT = visc.DiffVisc_dT(planet->phi(J-1) , planet->nu(J-1)*planet->rho(J-1) , planet->ym(J-1,0) , planet->ym(J-1,2) , planet->diffphiT(J-1));

}

EoS::EoS(double *_variables):
/* Visc params: phi_c , alpha , visc_l , visc_s(0) [ Pa s] , (T0 , V0 [m^3mol-1], E0 [ Jmol-1]) , smoothing parameter , visc limit (for smoothing) , Arrhenius on or off*/
    visc({0.4 , 26, 100.0, 1e21, 1600, 5e-6 , 300e3 , 3000.0 , 1e6 , -1})
{
     EoStype = "Mix";
    
    std::string filebase , x_name , y_name , f_name;
    int grid_size[2];
    grid_size[0] = 300;
    grid_size[1] = 300;

    x_name = "P";
    y_name = "T";
    filebase = "eos/Stixrude";
   
    std::cout << "Preparing solid EoS...\n" << "From folder:\n" << filebase << std::endl;

    f_name = "300rho";
    std::cout << "Preparing density..." << std::endl;
    interp_rho_s.setup_from_coeff(filebase, x_name, y_name, f_name, grid_size);
    
    f_name = "300del_rho";
    std::cout << "Preparing delta..." << std::endl;
    interp_del_rho_s.setup_from_coeff(filebase, x_name, y_name, f_name, grid_size);
    
    f_name = "300Cp";
    std::cout << "Preparing heat capacity..." << std::endl;
    interp_Cp_s.setup_from_coeff(filebase, x_name, y_name, f_name, grid_size);

    filebase = "eos/RTpress";
    std::cout << "Preparing melt EoS...\n" << "From folder:\n" << filebase << std::endl;
    
    grid_size[0] = 250;
    grid_size[1] = 250;
    f_name = "250rho";
    std::cout << "Preparing density..." << std::endl;
    interp_rho_m.setup_from_coeff(filebase, x_name, y_name, f_name, grid_size);
    
    f_name = "250del_rho";
    std::cout << "Preparing delta..." << std::endl;
    interp_del_rho_m.setup_from_coeff(filebase, x_name, y_name, f_name, grid_size);
    
    f_name = "250Cp";
    std::cout << "Preparing heat capacity..." << std::endl;
    interp_Cp_m.setup_from_coeff(filebase, x_name, y_name, f_name, grid_size);

    if(_variables[0] >= 1){
        
        grid_size[0] = 300;
        grid_size[1] = 300;

        if( _variables[0] == 1){
            filebase = "eos/GammaIron";
        }else if( _variables[0] == 2){
            filebase = "eos/EpsilonIron";
        }
        std::cout << "Preparing solid core EoS...\n" << "From folder:\n" << filebase << std::endl;

        f_name = "300rho";
        std::cout << "Preparing density..." << std::endl;
        interp_rho_c.setup_from_coeff(filebase, x_name, y_name, f_name, grid_size);
    
        f_name = "300del_rho";
        std::cout << "Preparing delta..." << std::endl;
        interp_del_rho_c.setup_from_coeff(filebase, x_name, y_name, f_name, grid_size);
    
        f_name = "300Cp";
        std::cout << "Preparing heat capacity..." << std::endl;
        interp_Cp_c.setup_from_coeff(filebase, x_name, y_name, f_name, grid_size);
        
       
    }
    std::cout << "No EoS setup needed." << std::endl;
};

void TemporaryFunction ()
{
    Eigen::ArrayXd temp;
    double *a;
    EoS temp_EoS(a);
    P_P temp_P(0,0);
    temp_EoS.compute_values<Eigen::ArrayXd>(&temp_P);
    temp_EoS.compute_derivs<Eigen::ArrayXd>(&temp_P);
    
}