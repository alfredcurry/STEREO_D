#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include "timestep.hpp"
#include "Jacob_Maker.hpp"
#include "Matrix_reader.h"
#include "StructStruct.hpp"
#include "PhysicalConsts.h"
#include "EquationOfState.hpp"


int main(int arg , char** argv){
    //grid size , mass , EoS type , surface temp OR luminosity ,  opacity , "n" , # , mu_m , tmax , Pbound type , f_time
    int J = std::strtol(argv[1] , nullptr , 0);
    int I = 4;
    double acc = 1e-7 ;
    double tmax = std::strtol(argv[9] , nullptr , 0);
    std::string infolder , outfolder ; 

    outfolder =  "Results/time/IdealGas/";
    ArrayXXd mPr0;
    ArrayXd scale0; 
    struct P_P planet(I,J);
    stepper tStepper(I,J);

    //planet properties
    planet.M = std::strtod(argv[2] , nullptr ) * M_Earth;
    planet.M_core = 0;
    planet.R_core = 0;
    
    planet.n = std::strtod(argv[6] , nullptr );

    planet.nu = ArrayXd::Constant(J,5.0e11);

    //read in guesses
   
    char Tstring[100] , condstring[100] , nstring[100] , Resstring[100];
    int a;
    a = snprintf(Resstring, 100 , "Res_%.1d_" , J);
   
   

        planet.Tbound = 'b';
        planet.Pbound = *argv[10];
        planet.nabtype = 'a';
    
        infolder = "InitialGuesses/HotPolytrope";
        planet.mu_m = std::strtod(argv[8], nullptr );
        a = snprintf(nstring, 100 , "n%.1f_%.3g_" , planet.n , boltz_R/planet.mu_m) ;

        if(planet.n == 1.0){
            planet.conds[0] = 0;
            planet.conds[1] = 0;
            planet.conds[2] = 0;
    
            planet.Tsurf[0] = std::strtod(argv[4] , nullptr );
            planet.Psurf[0] = std::strtod(argv[5] , nullptr );
            a = snprintf(Tstring, 100 , "Ts%.2g_Ps%.2g" , 10.0, 1e6);
            a = snprintf(condstring, 100 , "k_%.2g_al%.1f_bt%.1f_" , 0.0 , 0.0 , 0.0);
        }
        else
        if(planet.Pbound == 'o' || planet.Pbound == 'c'){
            
            planet.opac = std::strtod(argv[5] , nullptr );
            a = snprintf(Tstring, 100 , "L%.2g_opac%.2g" , std::strtod(argv[4] , nullptr ), planet.opac);
            a = snprintf(condstring, 100 , "AD" );

        }else
        if(planet.Pbound == 'f'){
            planet.Tsurf[0] = std::strtod(argv[4] , nullptr );
            planet.Psurf[0] = std::strtod(argv[5] , nullptr );
            a = snprintf(Tstring, 100 , "Ts%.2g_Ps%.2g" , planet.Tsurf[0] , planet.Psurf[0]); 
            planet.conds[0] = 0;
            planet.conds[1] = 0;
            planet.conds[2] = 0;
            a = snprintf(condstring, 100 , "k_%.2g_al%.1f_bt%.1f_" , planet.conds[0] , planet.conds[1] , planet.conds[2]);
            
        }
       
    

    double EoSvariables[2];
    EoSvariables[0] = planet.mu_m ;
    EoSvariables[1] = planet.n ;
    EoS EoSIdeal(EoSvariables);
    outfolder.append(EoSIdeal.EoStype);

    std::string txtstring = ("./"+infolder+"/results1"+nstring+Resstring+condstring+Tstring+".txt");
    std::string scalestring = ("./"+infolder+"/scale1"+nstring+Resstring+condstring+Tstring+".txt");
    std::cout << txtstring << "\n" << scalestring << std::endl;
    char *txtfile = &txtstring[0];
    char *scalefile = &scalestring[0];
    
    mPr0 = readArray(txtfile);

    scale0 = readArray(scalefile);
    
    planet.y = mPr0.block(0,1,J,I);

    planet.scale = scale0.tail(I);
    planet.H = planet.scale(3)/planet.M;
    //m  grid 
    planet.m = mPr0.col(0);
    JacobMaker::prepm( &planet );
    if(planet.n == 1.0){
        planet.scale(3) = 1e12 ; //4*M_PI*planet.scale(1)*planet.scale(1)*Sig_Stefan*pow(planet.Tsurf[0],4.0);
        
        planet.H = planet.scale(3)/planet.M;
        planet.opac = 1.7e-6 ; //G_Newt * planet.M * 2.0/3 /( pow(planet.scale(1),2) * planet.Psurf[0]) ; 
        std::cout << planet.scale(3) << " " << planet.opac << std::endl;
    }
    if (planet.M != scale0(0)){
        std::cout << "Wrong mass of planet" << std::endl;
    }
    if(planet.Pbound == 'c'){
        planet.opac = 2.0/3*G_Newt * planet.M * planet.scale(1)*planet.scale(1)/planet.opac;
    }
    double f_time = std::strtod(argv[11] , nullptr );
    char log_or_lin = '2';
    
    tStepper.timestepper_fullH(&planet , &EoSIdeal ,tmax , 1e40 ,0 , acc , f_time , 10 , argv[7] , &log_or_lin , outfolder );
    std::cout << "Completed" << std::endl;
    return 0;
}
