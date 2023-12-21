#define BOOST_TEST_MODULE Struct_test_boost
#include <boost/test/included/unit_test.hpp>
#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include "timestep.hpp"
#include "Jacob_Maker.hpp"
#include "Matrix_reader.h"
#include "StructStruct.hpp"
#include "PhysicalConsts.h"
#include "Interpolator.hpp"
#include "EquationOfState.hpp"

double Time_Struct(int _J , double n , int _tmax, double _ftime){
    int J = _J;
    int I = 4;
    double acc = 1e-7 ;
    double tmax = _tmax;
    std::string infolder , outfolder ; 

    outfolder =  "Tests/RunFiles/";
    ArrayXXd mPr0;
    ArrayXd scale0; 
    struct P_P planet(I,J);
    stepper tStepper(I,J);

    //planet properties
    planet.M = M_Earth;
    planet.M_core = 0;
    planet.R_core = 0;
    planet.L_core = 0;
    planet.n = n;
    planet.mu_m = 1.0;
    planet.nu = ArrayXd::Constant(J,5.0e11);
    
    planet.BL = '0';
    planet.MassLossType = '0';
    planet.RadioType = '0';
    double EoSvariables[2];
    EoSvariables[0] = planet.mu_m ;
    EoSvariables[1] = planet.n ;
    EoS EoSIdeal(EoSvariables);
    outfolder.append(EoSIdeal.EoStype);

    //read in guesses
   
    char Tstring[100] , condstring[100] , nstring[100] , Resstring[100];
    int a;
    
    planet.Tbound = 'b';
    planet.Pbound = 'c';
    planet.nabtype = 'a';
    
    infolder = "InitialGuesses/HotPolytrope";
    
    a = snprintf(Resstring, 100 , "Res_%.1d_" , J);
    a = snprintf(nstring, 100 , "n%.1f_%.3g_" , planet.n , boltz_R/planet.mu_m) ;

    
    planet.conds[0] = 0;
    planet.conds[1] = 0;
    planet.conds[2] = 0;
    
    a = snprintf(Tstring, 100 , "Ts%.2g_Ps%.2g" , 10.0, 1e6);
    a = snprintf(condstring, 100 , "k_%.2g_al%.1f_bt%.1f_" , 0.0 , 0.0 , 0.0);
     
    std::string txtstring = ("./"+infolder+"/results1"+nstring+Resstring+condstring+Tstring+".txt");
    std::string scalestring = ("./"+infolder+"/scale1"+nstring+Resstring+condstring+Tstring+".txt");

    std::cout << txtstring << "\n" << scalestring << std::endl;
    char *txtfile = &txtstring[0];
    char *scalefile = &scalestring[0];

    
    mPr0 = readArray(txtfile);
    std::cout << "will read scale0" << std::endl;

    scale0 = readArray(scalefile);
    std::cout << "scale0" << std::endl;

    std::cout << scale0 << std::endl;
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
    double f_time = _ftime;
    char num = 'T';
    char log_or_lin = '2';

    tStepper.timestepper_fullH(&planet , &EoSIdeal , tmax , 1e20 , 0 , acc , f_time , 11 , &num , &log_or_lin , outfolder );
    std::cout << "Model run. Checking against theoretical" << std::endl;
    //check against theoretical    
    ArrayXXd bulk_calc;
    double B , Ts0 , R0 , alpha = 4;
    ArrayXd t, Ks , R_th;

    std::cout << tStepper.print.bulkfile[0] << std::endl;
    char *bulkchar = &tStepper.print.bulkfile[0];
    bulk_calc = readArray(bulkchar);
    std::cout << "array read" << std::endl;
    R0 = bulk_calc(0,2);
    Ts0 = bulk_calc(0,6);
    t = bulk_calc.col(0) *yrs_to_s;
    B = G_Newt*M_Earth*M_Earth/2 - M_Earth*boltz_R*Ts0*R0 ;
    std::cout << B << std::endl;
    
    //Ks = pow(bulk_calc.col(6)*boltz_R,planet.n+1)/pow(bulk_calc.col(5),1/planet.n);
    R_th = R0 - 4*M_PI*Sig_Stefan*pow(Ts0*R0,4)*(t-t(0))/B ;
    std::cout << "rth " <<  R_th<< std::endl;
    //bulk_th.col(2) = bulk_calc(0,3)*R0/bulk_th.col(1);
   // bulk_th.col(0) = bulk_calc(0,1)*pow(R0/bulk_th.col(1),4);
    //bulk_th.col(5) = bulk_calc(0,6)*R0/bulk_th.col(1);
    //bulk_th.col(3) = 4*np.pi*bulk_th.col(1)*bulk_th.col(1)*pow(bulk_th.col(5),4);
    //bulk_th.col(4) = bulk_calc(0,5)*pow(R0/bulk_th.col(1),alpha);
    //bulk_th.col(6) = Ks(0)*R0**(1-3.0/planet.n)*pow(bulk_th.col(1),(-1+3/planet.n));

    std::cout << "Time_Poly_Test completed\n" << std::endl;

    return sqrt(pow((R_th-bulk_calc.col(2))/R_th,2.0).sum()/bulk_calc.rows());
}

BOOST_AUTO_TEST_CASE( Time_Poly_Test)
{
  boost::unit_test::unit_test_log.set_stream( std::cerr );
  boost::unit_test::unit_test_log.set_threshold_level(boost::unit_test::log_level::log_test_units); 
  
  BOOST_TEST(Time_Struct(500 , 1 , 50 , 0.01) < 3.2e-2);
  //for reference
  //[0.00033 , 0.001 , 0.0033 , 0.01 , 0.033 ]
  //[0.00097923 0.00328918 0.01135968 0.03155371 0.13581419]
}