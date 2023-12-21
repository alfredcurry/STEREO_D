#include "BoundaryCondition.hpp"
#include <vector>
#include "PhysicalConsts.h"

template<class Real>
Real B_C::T_surf(struct P_P *planet , Real &L){

    if(planet->Tbound != 't'){
        std::cout << "Tbound " << planet->Tbound << std::endl;
        throw std::invalid_argument("\nTemperature Boundary Condition Type Not Tabulated" );
    }
    int *indexes;
    Real L_tild = L/(2*M_PI*planet->scale(1)*planet->scale(1));
    //Real g = G_Newt*M/(R*R);
    std::cout << L <<" LR "<< planet->scale(1) <<" L tilda " << L_tild << " g "<< planet->g_surf << std::endl;
    indexes = interp_T0.return_index(log10(L_tild),planet->g_surf);  
    return interp_T0.return_value(log10(L_tild) , planet->g_surf , indexes[0] , indexes[1]);
}
template<class Real>
void B_C::T_surf_derivs(struct P_P *planet ){
    
    int *indexes;
    double *grads;
    Real L_tild = planet->scale(3)/(2*M_PI*planet->scale(1)*planet->scale(1));
    //Real g = G_Newt*M/(R*R);
    
    indexes = interp_T0.return_index(log10(L_tild),planet->g_surf);  
    grads = interp_T0.return_grads(log10(L_tild) , planet->g_surf , indexes[0] , indexes[1]);
    planet->Tsurf[1+1*3] = grads[0]/(log(10)*L_tild)/(2*M_PI*planet->scale(1)*planet->scale(1));
    planet->Tsurf[1+1*1] = grads[0]/(log(10)*L_tild)*-planet->scale(3)/(M_PI*pow(planet->scale(1),3)) ;//+ grads[1]*-2*G_Newt*M/(R*R*R);
    planet->Tsurf[1+1*0] = 0 ;
}

void B_C::read_table(double &P0 , double &Tss){}

B_C::B_C(){

    std::string filebase , x_name , y_name , f_name;
    int grid_size[2];
    grid_size[0] = 1000;
    grid_size[1] = 500;

    x_name = "logL";
    y_name = "g";
    filebase = "Boundary/TPvisc";//No_redist_P01e9_limited_Tss2145"; //
   
    std::cout << "Preparing Boundary condition...\n" << "From folder:\n" << filebase << std::endl;

    f_name = "T0";
    std::cout << "Preparing T0..." << std::endl;
    interp_T0.setup_from_coeff(filebase, x_name, y_name, f_name, grid_size);
}

void TemporaryFunction_BC()
{
    P_P planet_temp(1,1);
    double L;
    B_C temp_BC;
    temp_BC.T_surf<double>(&planet_temp,L);
    temp_BC.T_surf_derivs<double>(&planet_temp);

}