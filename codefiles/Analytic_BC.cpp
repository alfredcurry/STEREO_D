#include "BoundaryCondition.hpp"
#include <vector>
#include "PhysicalConsts.h"
#include "SurfaceTemp.hpp"
#include "Vector_reader.h"

/// Implementation of the B_C class (BoundaryConditon) for analytical formulations

template<class Real>
Real B_C::T_surf(struct P_P *planet , Real &L){

    if( planet->Tbound == 'f' ){  
        return planet->Tsurf[0];  
    }else{
        int J = planet->rho.size();
        if( planet->Tbound == 'b' ){
            return T_surface<double>::BB_temp(planet->scale(3)*planet->y(J-1,3) , planet->scale(1)*planet->y(J-1,1));
        }else if( planet->Tbound == 'i'){
            planet->Teff = T_surface<double>::BB_temp(planet->scale(3)*planet->y(J-1,3) , planet->scale(1)*planet->y(J-1,1));
            return T_surface<double>::IRR_temp(planet->Teff , planet->Teq , planet->tau , 1.0);
        }else if( planet->Tbound == 'm' ){
            planet->Teff = T_surface<double>::BB_temp(planet->scale(3)*planet->y(J-1,3) , planet->scale(1)*planet->y(J-1,1));
            return T_surface<double>::IRR_temp(planet->Teff , planet->Teq , 1.0/3*planet->tau , 0.5);
        }else if( planet->Tbound == 'l' ){
            int J = planet->rho.size();
            planet->Tsurf[0] = fit(planet->ym(J-1,3)/(2*M_PI*planet->ym(J-1,1)*planet->ym(J-1,1)) , planet->g(J-1));
            
            return planet->Tsurf[0];
	    }else {
            std::cout << "Tbound " << planet->Tbound << std::endl;
            throw std::invalid_argument("\nTemperature Boundary Condition Type Invalid" );
        }
    }
}

template<class Real>
void B_C::T_surf_derivs(struct P_P *planet){
    if( planet->Tbound == 'b' ){
        planet->Tsurf[1+1*1] = T_surface<double>::BB_diff_R(planet->scale(1),planet->Tsurf[0]);
        planet->Tsurf[1+1*3] = T_surface<double>::BB_diff_L(planet->scale(3),planet->Tsurf[0]);
    } else if (planet->Tbound == 'i' ){
        planet->Tsurf[1+1*1] = T_surface<double>::IRR_diff_R(planet->scale(1), planet->Tsurf[0] , planet->Teq , planet->Teff , planet->tau , 1.0);
        planet->Tsurf[1+1*3] = T_surface<double>::IRR_diff_L(planet->scale(3), planet->Tsurf[0] , planet->Teq , planet->Teff , planet->tau , 1.0);
    }else if (planet->Tbound == 'm' ){
        planet->Tsurf[1+1*1] = T_surface<double>::IRR_diff_R(planet->scale(1), planet->Tsurf[0] , planet->Teq , planet->Teff , planet->tau*1.0/3 , 0.5);
        planet->Tsurf[1+1*3] = T_surface<double>::IRR_diff_L(planet->scale(3), planet->Tsurf[0] , planet->Teq , planet->Teff , planet->tau*1.0/3 , 0.5);
    }else if( planet->Tbound == 'f' ){
        planet->Tsurf[1+1*1] = 0;
        planet->Tsurf[1+1*3] = 0;
    }else if( planet->Tbound == 'l' ){
        double T_dg , T_dF;
        T_dg = fit.g_derivative(planet->scale(3)/(2*M_PI*planet->scale(1)*planet->scale(1)) , planet->g(planet->g.size()-1));
        T_dF = fit.F_derivative(planet->scale(3)/(2*M_PI*planet->scale(1)*planet->scale(1)) , planet->g(planet->g.size()-1));
        //std::cout << planet->scale(3)/(2*M_PI*planet->scale(1)*planet->scale(1)) << " " << planet->g(planet->g.size()-1) << std::endl;
        //std::cout << "surf dervis" << T_dg << " " << T_dF << std::endl;
        planet->Tsurf[1+1*0] = 0;
        planet->Tsurf[1+1*1] = -(T_dg * 2.0*G_Newt*planet->M + T_dF * planet->scale(3)/M_PI )/pow(planet->scale(1),3);
        planet->Tsurf[1+1*2] = 0;
        planet->Tsurf[1+1*3] = T_dF / (2*M_PI*planet->scale(1)*planet->scale(1));
    }else {
        std::cout << "Tbound " << planet->Tbound << std::endl;
        throw std::invalid_argument("\nTemperature Boundary Condition has no derivative" );
        }
}

/// @brief Read in parameters for the fitted function for the tidally locked BC
/// @param P0 pressure at edge of 1D grid
/// @param Tss substellar temperature
void B_C::read_table(double &P0,double &Tss){
    //double A1 =  8.85558560e+04 , a1 = 1.0 , b1 = -1.0  , A2 = 1.69223668e+03 , a2 = 2.59215776e+01 , b2 = 1.09102064e+01 , A4 = 0.5 , a4 = 0.5 , b4 = -2.44372441e-01 , beta = 8.33435231e+00 , alpha = 1.88064746e+01 , Tcr , gam = 12.0;
    std::vector<double> variables(14);
    std::string parameter_file ;
   // variables = { A1 , a1 , b1 , A2 , a2 , b2 , A4 , a4 , b4 , beta , alpha , gam , Tcr };
    
    bool low_mass = true, inward = false , spherical = false;
    if(low_mass == true && inward == false && spherical == false ){
        char Tssstring[20];
        parameter_file = "./Boundary/fitted_parameters/";
        
        snprintf(Tssstring, 20 , "P0%.0g/Tss%.0f.txt" ,P0 , Tss);
        
        parameter_file = parameter_file+Tssstring;
        std::cout << "Read fitted parameters" << std::endl;
        readVector(variables , 14 , parameter_file);
        
    }else
    if(low_mass == true && inward == true && spherical == false ){
        char Tssstring[20];
        parameter_file = "./Boundary/fitted_parameters_inward/";
        snprintf(Tssstring, 20 , "P0%.0g/Tss%.0f.txt" ,P0 , Tss);
        
        parameter_file = parameter_file+Tssstring;
        std::cout << "Read fitted parameters (inward included)" << std::endl;
        readVector(variables , 14 , parameter_file);

    }else
    if(low_mass == true && inward == false && spherical == true ){
        char Tssstring[30];
        parameter_file = "./Boundary/fitted_parameters/";
        snprintf(Tssstring, 30 , "P0%.0g/spherical_Teq%.0f.txt" ,P0 , Tss);
        
        parameter_file = parameter_file+Tssstring;
        std::cout << "Read fitted parameters (spherical)" << std::endl;
        readVector(variables , 14 , parameter_file);

    }else{
        parameter_file = parameter_file+"Teq"+std::to_string(300)+".txt" ;
        readVector(variables , 13 , parameter_file);
    }
    

    std::cout <<"varibales " <<  variables.size() << std::endl;
    for(int i=0 ; i < variables.size() ; i++){
        std::cout << variables[i] << std::endl;
    }
    std::cout << "vars in" << std::endl;
    fit.Assign(variables);
}

B_C::B_C(){
}

void TemporaryFunction_BC()
{
    P_P planet_temp(1,1);
    double L;
    B_C temp_BC;
    temp_BC.T_surf<double>(&planet_temp,L);
    temp_BC.T_surf_derivs<double>(&planet_temp);

}