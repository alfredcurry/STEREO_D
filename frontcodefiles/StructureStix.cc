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
#include <stdexcept>

using namespace Eigen;

int main(int arg , char** argv){
    ///Run-time model parameters:
    ///grid size , mass , Tbound type , Pbound type , surface temp OR L, Psurf /opacity , # , ftime , tmax , folder , nabtype , MLtype , MassLoss , BLtype , corefrac , Tss
    
    ///Read in or assign parameters
    int J = std::strtol(argv[1] , nullptr , 0);
    int I = 4;
    double acc = 1e-6;
    std::string folder = "Results/" ,  eos_name(argv[11]);
    folder = folder+eos_name ; 
    struct P_P planet(I,J);
    //std::cout << "Jcore " << planet.Jcore << std::endl;
    planet.M = std::strtod(argv[2] , nullptr ) * M_Earth;
    planet.M_core = 0;
    planet.R_core = 0;
    planet.nabtype = *argv[12];
    planet.Teq = std::strtod(argv[17] , nullptr );
    planet.tau = 0;
    planet.a_crystal = 1e-3; //crystal size / m
    planet.BL = *argv[15];
    planet.RadioType = 'E';

    planet.Tbound = *argv[3];
    planet.Pbound = *argv[4];
    planet.ML = *argv[13];
    planet.MassLossType = *argv[14] ; 
    
    double EoSvariables[1];

    planet.P_min = std::strtod(argv[6] , nullptr ); 
    planet.Psurf[0] = planet.P_min;
    planet.Tsurf[0] = std::strtod(argv[5] , nullptr );
    
    planet.conds[0] = 4.0; // Bower (2017) 4.0 W m-1 K-1
    planet.conds[1] = 0;
    planet.conds[2] = 0;

    double core_frac = std::strtod(argv[16] , nullptr ) , Rcore ;
    
    if(core_frac == 0){
        EoSvariables[0] = 0;
    }else if ( planet.M <= 0.3* M_Earth){
        EoSvariables[0] = 1;
    }else{
        EoSvariables[0] = 2;
    }

    EoS EoSStix(EoSvariables);
    char* num = argv[7];
    char *txtfile , *scalefile;

    ArrayXXd Allthing(J , I + I + 3) , Fls(J,3);
    
    stepper tStepper(I,J);
    if(planet.Tbound == 'l'){
        tStepper.con.Jacob.boundary.read_table(planet.Psurf[0],planet.Teq);
    }

    planet.m = JacobMaker::power_law_m(J,1,1.3);//  ArrayXd::LinSpaced(J,0,1);//::massbins(J);
    
    //set up intitial guess of structure (read from file)
    ArrayXXd mPr0;
    ArrayXd scale0; 
    

    JacobMaker::prepm( &planet );
  
    std::cout << J << " " << planet.M << " " << planet.Tsurf[0] <<" " << planet.n << std::endl;  
    if(planet.M <= 1.1*M_Earth && planet.M >0.55*M_Earth){    
        std::cout << "Mass " << planet.M /M_Earth << std::endl; 
        mPr0 = readArray("./InitialGuesses/MgSiO3perov/results1MgSiO3_M_E1_Res_1000_ADTs3e+02_Ps1e+05.txt");
        scale0 = readArray("./InitialGuesses/MgSiO3perov/scale1MgSiO3_M_E1_Res_1000_ADTs3e+02_Ps1e+05.txt");
    }else if(planet.M > 0.075*M_Earth && planet.M < 0.15*M_Earth){
        std::cout << "Mass " << planet.M /M_Earth << std::endl; 
        mPr0 = readArray("./InitialGuesses/MgSiO3perov/results1MgSiO3_M_E0.1_Res_600ADTeq2000tau0_Ps1.txt");
        scale0 = readArray("./InitialGuesses/MgSiO3perov/scale1MgSiO3_M_E0.1_Res_600ADTeq2000tau0_Ps1.txt");
        //mPr0 = readArray("./InitialGuesses/MgSiO3perov/results_transition.txt");
        //scale0 = readArray("./InitialGuesses/MgSiO3perov/scale_transition.txt");
    }else if(planet.M <= 0.25*M_Earth && planet.M >= 0.15*M_Earth ){
        std::cout << "Mass " << planet.M /M_Earth << std::endl; 
        mPr0 = readArray("./InitialGuesses/MgSiO3perov/results1MgSiO3_E_M_E0.2_Res_800_ADL0_Teq2145tau0Ps1e+03.txt");
        scale0 = readArray("./InitialGuesses/MgSiO3perov/scale1MgSiO3_E_M_E0.2_Res_800_ADL0_Teq2145tau0Ps1e+03.txt");
        
    }else if(planet.M <= 0.55*M_Earth && planet.M > 0.25*M_Earth ){
        std::cout << "Mass " << planet.M /M_Earth << std::endl; 
        mPr0 = readArray("./InitialGuesses/MgSiO3perov/results1MgSiO3_E_M_E0.3_Res_800_ADL0_Teq2145tau0Ps1e+03.txt");
        scale0 = readArray("./InitialGuesses/MgSiO3perov/scale1MgSiO3_E_M_E0.3_Res_800_ADL0_Teq2145tau0Ps1e+03.txt");

    }else if(planet.M == 0.03*M_Earth || planet.M == 0.035*M_Earth || planet.M == 0.04*M_Earth){
        std::cout << "Mass " << planet.M /M_Earth << std::endl; 
        mPr0 = readArray("./InitialGuesses/MgSiO3perov/results1MgSiO3_E_M_E0.03_Res_400_ADTs2200Ps1e+08.txt");
        scale0 = readArray("./InitialGuesses/MgSiO3perov/scale1MgSiO3_E_M_E0.03_Res_400_ADTs2200Ps1e+08.txt");
    }else if(planet.M >= 0.04*M_Earth && planet.M <= 0.075*M_Earth){
        std::cout << "Mass " << planet.M /M_Earth << std::endl; 
        mPr0 = readArray("./InitialGuesses/MgSiO3perov/results1MgSiO3_E_M_E0.05_Res_400_ADTs2000Ps1e+09.txt");
        scale0 = readArray("./InitialGuesses/MgSiO3perov/scale1MgSiO3_E_M_E0.05_Res_400_ADTs2000Ps1e+09.txt");
    }else{
        throw std::invalid_argument("Not a mass with a guess");
    }

    ArrayXd mm0 = JacobMaker::diffminus(mPr0.col(0));
    double mtemp = planet.mm(0);
    planet.mm(0) = 0;

    mm0 = mPr0.col(0) - mm0/2 ;
    mm0(0) = 0;
    mm0(mm0.size()-1) = 1;

    planet.y.col(0) = interpolate(planet.mm , mm0 , mPr0.col(1) ); // pressure, using interpolation
    planet.y.block(0,1,J-1,1) = interpolate(planet.m.tail(J-1) , mPr0.col(0) , mPr0.col(2) ); //radius , using interpolation
    planet.mm(0) = mtemp;
 
    if( abs(planet.y(J-3,1) - planet.y(J-2,1))/planet.y(J-2,1) < 1e-6){
        planet.y(J-3,1) = (planet.y(J-4,1) + planet.y(J-2,1))/2;
    }
    planet.y.block(0,3,J-1,1) = planet.m.tail(J-1); 

    double Tc0 , rhoc0;
    int a;
    
    planet.scale(0) = scale0(1);
    planet.scale(1) = scale0(2);
    planet.scale(2) = planet.Tsurf[0] + 1e-7*pow(planet.Tsurf[0],3) ;    
    planet.y.col(2) = ArrayXd::Constant(J , 1.0 ) ;
    planet.scale(3) = 4*M_PI * planet.scale(1)*planet.scale(1) * std::strtod(argv[10] , nullptr ); 
    planet.H = planet.scale(3)/planet.M;
    
    planet.y(J-1,1) = planet.y(J-2,1);
    planet.y(J-1,3) = planet.y(J-2,3);
    planet.y(J-1,0) = planet.Psurf[0]/planet.scale(0);

    double kaim = planet.conds[0];
    std::cout << "\n" << planet.scale << std::endl;
    std::cout << "Tsurf " << planet.Tsurf[0] << " T " << planet.y(J-1,2)* planet.scale(2) << " P " << planet.Psurf[0] << std::endl;

    planet.T0 = planet.scale(2)*planet.y.col(2);
    planet.P0 = planet.scale(0)*planet.y.col(0);

    //Prepare the timestepper and run it!
    int t = 0;
    double M_step = planet.M*0.01 , Maim = 10.0*planet.M;
    planet.dt = 1e308;

    int tmax = 100000;
    double timemax = std::strtod(argv[9] , nullptr );

    double ftime = std::strtod(argv[8] , nullptr );
    double Taim = planet.Tsurf[0] ; //std::strtod(argv[8] , nullptr );

    Rcore = pow(0.75 / M_PI * core_frac * planet.M / 13000 , 1.0/3);
    planet.L_core = 0;//planet.scale(3);

    char log_or_lin = '2';
    //tStepper.core_incrementer2(&planet , &EoSStix , Rcore , core_frac, acc , num , folder);
    planet.g_surf = G_Newt*planet.M/(planet.scale(1)*planet.scale(1));

    std::cout << "About to start" << std::endl;
    tStepper.introduce_core(&planet , &EoSStix , core_frac, acc , num  , folder);
    //tStepper.L_alter(&planet , &EoSStix , 2e15 , acc , num  , folder);
    //tStepper.temp_alter(&planet , &EoSStix , Taim , acc , num , folder);
    planet.g_surf = G_Newt*planet.M/(planet.scale(1)*planet.scale(1));
    //tStepper.simp_M_alter(&planet , &EoSStix , Maim , 0.1 , acc , num , folder);
    
    if( planet.MassLossType != '0' && ftime == 1){
        tStepper.static_mass_loss(&planet , &EoSStix , tmax , timemax , acc , 1 , num , &log_or_lin , folder);
    }else{
        tStepper.timestepper_fullH(&planet , &EoSStix , tmax , timemax , 0 , acc , ftime , 1 , num , &log_or_lin , folder);
    }
    //tStepper.timestepper_Hdecay(&planet , &EoSStix , tmax , 1e10 , 1e13 , acc , ftime , 10 , num , folder);

    return 0;
}
