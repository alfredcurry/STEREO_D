#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include "Matrix_reader.h"
#include <string>
#include <vector>
#include "Interpolator.hpp"
#include "PhysicalConsts.h"
#include "Jacob_Maker.hpp"
#include "SurfaceTemp.hpp"
#include "SurfacePressure.hpp"
#include "EquationOfState.hpp"
#include "timestep.hpp"

using namespace Eigen;

int main(int arg , char** argv){
    //grid size , mass , Tbound type , Pbound type , surface temp OR L, Psurf /opacity , # , ftime , tmax , folder
    int J = std::strtol(argv[1] , nullptr , 0);
    int I = 4;
    double acc = 1e-5;
    std::string folder = "Results/" ,  eos_name(argv[10]);
    folder = folder+eos_name ; 
    char Tstring[100] , condstring[100] , Resstring[100], massstring[100];
    struct P_P planet(I,J);
    planet.M = std::strtod(argv[2] , nullptr ) * M_Earth;
    
    double core_frac = 0.326 , Rcore = 3.5e6;
    planet.M_core = core_frac*planet.M;
    planet.R_core = Rcore;
    planet.nabtype = 'a';
    planet.Teq = 1000;

    planet.Tbound = *argv[3];
    planet.Pbound = *argv[4]; 

    double EoSvariables[1];
    EoSvariables[0] = 0;
   
    EoS EoSStix(EoSvariables);
    char* num = argv[7];
    char *txtfile , *scalefile;

    ArrayXXd Allthing(J , I + I + 3) , Fls(J,3);

    //m and guesses

    planet.m = JacobMaker::power_law_m(J , (planet.M - planet.M_core)/planet.M , 1.0 );
  
    stepper tStepper(I,J);
    ArrayXXd mPr0;
    ArrayXd scale0; 
    //set up guess
    JacobMaker::prepm( &planet );
  
    //std::cout << J << " " << planet.M <<" " << planet.EoS << " " << planet.Tsurf[0] <<" " << planet.n << std::endl;  
        
  //  mPr0 = readArray("./InitialGuesses/MgSiO3perov/results1MgSiO3_M_E1_Res_1000_ADTs3e+02_Ps1e+05.txt");
//    scale0 = readArray("./InitialGuesses/MgSiO3perov/scale1MgSiO3_M_E1_Res_1000_ADTs3e+02_Ps1e+05.txt");

    mPr0 = readArray("./InitialGuesses/MgSiO3perov/results1MgSiO3_M_E1_Res_1000Mc1.95e+24_Rc3.5e+06_Lc1.58e+20ADTs1500Ps1e+09.txt");
    scale0 = readArray("./InitialGuesses/MgSiO3perov/scale1MgSiO3_M_E1_Res_1000Mc1.95e+24_Rc3.5e+06_Lc1.58e+20ADTs1500Ps1e+09.txt");


    ArrayXd mm0 = JacobMaker::diffminus(mPr0.col(0));
    double mtemp = planet.mm(0);
    planet.mm(0) = 0;

    mm0 = mPr0.col(0) - mm0/2 ;
    mm0(0) = 0;
    mm0(mm0.size()-1) = 1;

    planet.y.col(0) = interpolate(planet.mm , mm0 , mPr0.col(1) ); // pressure, using interpolation
    planet.y.block(0,1,J-1,1) = interpolate(planet.m.tail(J-1) , mPr0.col(0) , mPr0.col(2) ); //radius , using interpolation
    planet.mm(0) = mtemp;
 
    planet.y.block(0,3,J-1,1) = ArrayXd::Constant(J-1,1.0); 

    double Tc0 , rhoc0;
    int a;

    
    planet.y.col(2) = ArrayXd::Constant(J,1.0);
    planet.Psurf[0] = std::strtod(argv[6] , nullptr ); 
    planet.Tsurf[0] = std::strtod(argv[5] , nullptr );
    planet.scale(0) = scale0(1);
    planet.scale(1) = scale0(2);
    planet.scale(2) = planet.Tsurf[0];    
    planet.scale(3) = 4*M_PI * planet.scale(1)*planet.scale(1) *Sig_Stefan * (pow(planet.Tsurf[0],4) - pow(planet.Teq,4));//(planet.H*planet.dm).sum() * planet.M;
    planet.H = 0; //planet.scale(3)/planet.M;

    planet.conds[0] = 0;
    planet.conds[1] = 0;
    planet.conds[2] = 0;
            
    a = snprintf(Tstring, 100 , "Ts%.2g_Ps%.2g" , planet.Tsurf[0] , planet.Psurf[0]);  
    
    if(planet.nabtype == 'a'){
        a = snprintf(condstring, 100 , "AD" );
    }else if (planet.nabtype == 'i'){
        a = snprintf(condstring, 100 , "ISO" );
    }else{    
        a = snprintf(condstring, 100 , "k_%.2g_al%.1f_bt%.1f_" , planet.conds[0] , planet.conds[1] , planet.conds[2]);
    }
    planet.y(J-1,1) = planet.y(J-2,1);
    planet.y(J-1,3) = planet.y(J-2,3);
    planet.y(J-1,0) = planet.Psurf[0]/planet.scale(0);
    planet.y(J-1,2) = planet.Tsurf[0]/planet.scale(2);
    //1.01*pow(planet.scale(3)/(4.0*M_PI*Sig_Stefan*planet.scale(1)*planet.scale(1)),1.0/4)/planet.scale(2);

    double kaim = planet.conds[0];
    std::cout << planet.y << std::endl;
    std::cout << "\n" << planet.scale << std::endl;
    std::cout << planet.Tsurf[0] << " T P " << planet.Psurf[0] << std::endl;
    planet.nu = ArrayXd::Constant(J,5.0e11);

    planet.T0 = planet.scale(3)*planet.y.col(2);
    planet.P0 = planet.scale(0)*planet.y.col(0);


    int t = 0;
    double M_step = planet.M*0.01 , Maim = planet.M;
    planet.dt = 1e308;

    int tmax = std::strtol(argv[9] , nullptr , 0);
    double ftime = std::strtod(argv[8] , nullptr );
    double Taim = std::strtod(argv[8] , nullptr );

    planet.L_core = planet.scale(3);

    planet.M_core = core_frac*planet.M;
    planet.R_core = Rcore;
    tStepper.temp_alter(&planet , &EoSStix , Taim , acc , num , folder);
    //tStepper.L_core_cheat(&planet , &EoSStix , 1e20 , acc , num , folder);
    
    //tStepper.simp_M_alter(&planet , &EoSStix , Maim , 0.1 , acc , num , folder);
    //tStepper.timestepper_full(&planet , &EoSStix , tmax , 1e10 , 2e19 , acc , ftime , 10 , num , folder);
    //tStepper.timestepper_Hdecay(&planet , &EoSStix , tmax , 1e10 , 1e13 , acc , ftime , 10 , num , folder);

    return 0;
}
