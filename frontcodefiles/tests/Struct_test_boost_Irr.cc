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

double IrrEnv_Struct(int J , double ftime){
     //grid size , # , tmax , ftime , tau , Teq , Lcore 
    
    int I = 4;
    int a;
    double acc = 1e-5 ;
    int tmax = 1000;
    double yrs_max = 1e10; 
    std::string outfolder ; 
    char *txtfile , *scalefile;
    outfolder =  "Tests/RunFiles/";
    ArrayXXd mPr0;
    ArrayXd scale0; 
	ArrayXXd Allthing(J,I+1);    
	struct P_P planet(I,J);
    stepper tStepper(I,J);

    //planet properties
    planet.M = (1.0/0.9) *  M_Earth;
    //planet.EoS = 'I' ; //s - Seagar(2007) (for perovskite) ; I - ideal gas , R = 8300
    planet.n = 2.5;
    planet.mu_m = 2;
    
    planet.BL = '0';
    planet.MassLossType = '0';
    
    planet.M_core =  M_Earth;
    planet.R_core = pow(planet.M_core/M_Earth,0.25)*R_Earth;
	double Lcore_f = 0; 
    
    double EoSvariables[2];
    EoSvariables[0] = planet.mu_m ;
    EoSvariables[1] = planet.n ;
    EoS EoSIdeal(EoSvariables);
    outfolder.append(EoSIdeal.EoStype);

    double alpha_Owen , beta_Owen , kap0_Owen;
    alpha_Owen = 0.68;
    beta_Owen = 0.45;
    kap0_Owen = 2.29e-8; //in SI : cgs - 0.1Pa 4.79e-8*10^alpha?  an alternative 2e-9 
    planet.Tbound = 'i';
    planet.Pbound = 'o';
    planet.nabtype = 'c';
    planet.ML = 'c';
    
    planet.nu = ArrayXd::Constant(J,1.0e-35);

    planet.Psurf[0] = 400;
    planet.Tsurf[0] = 300;

    double kaim = 16*Sig_Stefan/(3 * kap0_Owen) * pow(boltz_R/planet.mu_m, -alpha_Owen);
    planet.conds[0] = kaim;    
    planet.conds[1] = -1-alpha_Owen;        
    planet.conds[2] = 3-alpha_Owen-beta_Owen;
    
    //prepare mass bins
    planet.m = JacobMaker::power_law_Plus(J , (planet.M - planet.M_core)/planet.M , 2.0 , 0.6); // (planet.n+1)/(3-planet.n)
    JacobMaker::prepm( &planet );

    //read in guesses
    std::string txtstring , scalestring;

        
    txtstring = ("InitialGuesses/CoredPolytrope/results1n2.5_2_Res_1000Mc5.97e+24_Rc6.37e+06_Lc1e+16k0_0.0458_al-1.68_bt1.87L1e+16_Ts300Ps2.5e+03.txt");//Results/time/IrrEnv4/results1n2.5_2_Res_500Mc5.97e+24_Rc6.37e+06_Lc4e+18k0_0.0458_al-1.68_bt1.87L4e+18_Teq300tau0.66.txt");//InitialGuesses/CoredPolytrope/results1n2.5_2_Res_500Mc5.97e+24_Rc6.37e+06_Lc1e+21k0_0.0458_al-1.68_bt1.87L1e+21_Teq300tau0.66.txt");//InitialGuesses/CoredPolytrope/results1n2.5_2_Res_200Mc5.97e+24_Rc6.37e+06_Lc1e+18k0_151_al-1.00_bt1.00L1e+18_Teq300tau100.txt");//Results/time/IrrEnv4/results1n2.5_2_Res_300Mc5.97e+24_Rc6.37e+06_Lc1e+21k0_0.0458_al-1.68_bt1.87L1e+21_Teq300tau2.txt");//  InitialGuesses/CoredPolytrope/results1n2.5_8.3e+03_Res_400Mc5.97e+24_Rc6.37e+06_Lc1e+20k0_0.00286_al-1.68_bt1.87L1e+20_Ts300Ps5e+03.txt");//"Results/time/IrrEnv3/results1n2.5_8.3e+03_Res_1000_ftime_7e-310Mc1.19e+25_Rc7.58e+06_Lc3e+19k0_0.00286_al-1.68_bt1.87L3e+19_Teq55tau1t=0.txt") ;//"+infolder+"/results1"+nstring+Resstring+corestring+condstring+Tstring+".txt");
    scalestring = ("InitialGuesses/CoredPolytrope/scale1n2.5_2_Res_1000Mc5.97e+24_Rc6.37e+06_Lc1e+16k0_0.0458_al-1.68_bt1.87L1e+16_Ts300Ps2.5e+03.txt");//Results/time/IrrEnv4/scale1n2.5_2_Res_500Mc5.97e+24_Rc6.37e+06_Lc4e+18k0_0.0458_al-1.68_bt1.87L4e+18_Teq300tau0.66.txt");//InitialGuesses/CoredPolytrope/scale1n2.5_2_Res_500Mc5.97e+24_Rc6.37e+06_Lc1e+21k0_0.0458_al-1.68_bt1.87L1e+21_Teq300tau0.66.txt");//Results/time/IrrEnv4/scale1n2.5_2_Res_300Mc5.97e+24_Rc6.37e+06_Lc1e+21k0_0.0458_al-1.68_bt1.87L1e+21_Teq300tau2.txt");//  InitialGuesses/CoredPolytrope/scale1n2.5_2_Res_200Mc5.97e+24_Rc6.37e+06_Lc1e+18k0_151_al-1.00_bt1.00L1e+18_Teq300tau100.txt");//InitialGuesses/CoredPolytrope/scale1n2.5_8.3e+03_Res_400Mc5.97e+24_Rc6.37e+06_Lc1e+20k0_0.00286_al-1.68_bt1.87L1e+20_Ts300Ps5e+03.txt");//"Results/time/IrrEnv3/scale1n2.5_8.3e+03_Res_1000_ftime_7e-310Mc1.19e+25_Rc7.58e+06_Lc3e+19k0_0.00286_al-1.68_bt1.87L3e+19_Teq55tau1t=0.txt"); // ("./"+infolder+"/scale1"+nstring+Resstring+corestring+condstring+Tstring+".txt");
    std::cout << txtstring << std::endl;
	 
    txtfile = &txtstring[0];
    scalefile = &scalestring[0];

    mPr0 = readArray(txtfile);
    scale0 = readArray(scalefile);   

    ArrayXd mm0, m0;
	m0 = mPr0.col(0) + planet.M_core/planet.M ;
    
	mm0 = JacobMaker::diffminus(m0) ; 
    double mtemp = planet.mm(0);
    planet.mm(0) = planet.M_core/planet.M;
    mm0 = m0 - mm0/2 ;
    mm0(0) = planet.mm(0);
    mm0(mm0.size()-1) = 1;

    planet.y.col(0) = interpolate(planet.mm , mm0 , mPr0.col(1) ); // pressure, using interpolation
    planet.y.block(0,1,J-1,1) = interpolate(planet.m.tail(J-1) , m0 , mPr0.col(2) ); //radius , using interpolation
    planet.y(J-2,1) = 1;
	planet.y(J-1,1) = planet.y(J-2,1);
    planet.y.col(2) = interpolate(planet.mm , mm0 , mPr0.col(3) ); // temp, using interpolation
    planet.y.block(0,3,J-1,1) = interpolate(planet.m.tail(J-1) , m0 , mPr0.col(4) ); //L, using interpolation
    planet.y(J-2,3) = 1;
	planet.y(J-1,3) = planet.y(J-2,3);

    planet.mm(0) = mtemp;
	
    planet.scale = scale0.segment(1,I);
	
    planet.L_core = planet.scale(3); // 1400*

    planet.H = 0*planet.scale(3)/planet.M;

    planet.Teq = 300 ; 
    planet.tau = 0.66 ;
    if (abs(planet.M - scale0(0)) > planet.M *0.0001){
        std::cout << planet.M << " " << scale0(0) << std::endl;
          throw std::invalid_argument("Wrong mass!" );
    }
    char num = 'T';
    char log_or_lin = '2';

    tStepper.timestepper_full(&planet, &EoSIdeal , tmax , yrs_max,  Lcore_f  ,  acc , ftime , 1 , &num , &log_or_lin , outfolder);

    char Richardfile[100];
    ArrayXXd Best_Array;
    ArrayXd Best_m , Best_rho , my_rho;
    snprintf(Richardfile, 100 , "./Tests/Best/Booth_planet_cool_%.1gyr.txt", yrs_max);
    std::cout << Richardfile << std::endl;
    Best_Array = readArray(Richardfile);
    Best_m = Best_Array.col(1)/M_Earth/1000;
    Best_rho = Best_Array.col(4)*1000;
    my_rho = interpolate(Best_m , planet.m -0.9 , planet.rho);
    std::cout << planet.m.transpose() -0.9<< std::endl;
    std::cout << Best_m.transpose() << std::endl;

    double rms = sqrt(pow((Best_rho - my_rho)/my_rho,2).sum()/Best_rho.size());
    std::cout << "Irr_Env_Test Completed\n" << rms << std::endl;
    return rms;
}

BOOST_AUTO_TEST_CASE( Irr_Env_Test)
{
  boost::unit_test::unit_test_log.set_stream( std::cerr );
  boost::unit_test::unit_test_log.set_threshold_level(boost::unit_test::log_level::log_test_units); 
  
  double err = IrrEnv_Struct(500 , 0.01); 
  std::cout << err << std::endl;
  BOOST_TEST( err < 0.13);
  //for reference
  // [200 , 500 , 1000]
  // [0.21279502 0.11268531 0.05954107]
}

