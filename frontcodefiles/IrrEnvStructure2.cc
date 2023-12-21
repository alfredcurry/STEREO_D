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

int main(int arg , char** argv){
    //grid size , # , tmax , ftime , tau , Teq , Lcore 
    int J = std::strtol(argv[1] , nullptr , 0);
    int I = 4;
    int a;
    double acc = 1e-5 ;
    double yrs_max = std::strtod(argv[3] , nullptr );
    std::string infolder , outfolder ; 
    char *txtfile , *scalefile;
	double ftime = std::strtod(argv[4] , nullptr );
    outfolder =  "Results/time/IrrEnv/";
    infolder = "InitialGuesses/CoredPolytrope";
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
	double Lcore_f = std::strtod(argv[7] , nullptr ); 
    
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
    planet.m = JacobMaker::power_law_Plus(J , (planet.M - planet.M_core)/planet.M , 3.4 , 0.7); // (planet.n+1)/(3-planet.n)
    JacobMaker::prepm( &planet );

    //read in guesses
    std::string txtstring , scalestring;

        //char Tstring[100] , condstring[100] , nstring[100] , Resstring[100] , corestring[100];
        //a = snprintf(nstring, 100 , "n%.1f_%.3g_" , planet.n , planet.boltzrel) ;
        //a = snprintf(Tstring, 100 , "L%.2g_Ts%.3gPs%.2g" , planet.scale(3) , planet.Tsurf[0] , planet.Psurf[0] ); 
    
        //a = snprintf(Tstring, 100 , "Ts%.2g_Ps%.2g_L%.2g" , planet.Tsurf[0] , planet.Psurf[0] , planet.scale(3)); 
        //a = snprintf(Tstring, 100 , "Ps%.2g_L%.2g" , planet.Psurf[0] , planet.scale(3)); 
        //a = snprintf(corestring, 100 , "Mc%.3g_Rc%.3g" , planet.M_core , planet.R_core);  
        //a = snprintf(condstring, 100 , "Walpk0_%.3g_al%.2f_bt%.2f" , planet.conds[0], planet.conds[1] , planet.conds[2]);
        //a = snprintf(Resstring, 100 , "Res_%.1d" , 500);
        txtstring = ("InitialGuesses/CoredPolytrope/results1n2.5_2_Res_1000Mc5.97e+24_Rc6.37e+06_Lc1e+16k0_0.0458_al-1.68_bt1.87L1e+16_Ts300Ps2.5e+03.txt");//Results/time/IrrEnv4/results1n2.5_2_Res_500Mc5.97e+24_Rc6.37e+06_Lc4e+18k0_0.0458_al-1.68_bt1.87L4e+18_Teq300tau0.66.txt");//InitialGuesses/CoredPolytrope/results1n2.5_2_Res_500Mc5.97e+24_Rc6.37e+06_Lc1e+21k0_0.0458_al-1.68_bt1.87L1e+21_Teq300tau0.66.txt");//InitialGuesses/CoredPolytrope/results1n2.5_2_Res_200Mc5.97e+24_Rc6.37e+06_Lc1e+18k0_151_al-1.00_bt1.00L1e+18_Teq300tau100.txt");//Results/time/IrrEnv4/results1n2.5_2_Res_300Mc5.97e+24_Rc6.37e+06_Lc1e+21k0_0.0458_al-1.68_bt1.87L1e+21_Teq300tau2.txt");//  InitialGuesses/CoredPolytrope/results1n2.5_8.3e+03_Res_400Mc5.97e+24_Rc6.37e+06_Lc1e+20k0_0.00286_al-1.68_bt1.87L1e+20_Ts300Ps5e+03.txt");//"Results/time/IrrEnv3/results1n2.5_8.3e+03_Res_1000_ftime_7e-310Mc1.19e+25_Rc7.58e+06_Lc3e+19k0_0.00286_al-1.68_bt1.87L3e+19_Teq55tau1t=0.txt") ;//"+infolder+"/results1"+nstring+Resstring+corestring+condstring+Tstring+".txt");
        scalestring = ("InitialGuesses/CoredPolytrope/scale1n2.5_2_Res_1000Mc5.97e+24_Rc6.37e+06_Lc1e+16k0_0.0458_al-1.68_bt1.87L1e+16_Ts300Ps2.5e+03.txt");//Results/time/IrrEnv4/scale1n2.5_2_Res_500Mc5.97e+24_Rc6.37e+06_Lc4e+18k0_0.0458_al-1.68_bt1.87L4e+18_Teq300tau0.66.txt");//InitialGuesses/CoredPolytrope/scale1n2.5_2_Res_500Mc5.97e+24_Rc6.37e+06_Lc1e+21k0_0.0458_al-1.68_bt1.87L1e+21_Teq300tau0.66.txt");//Results/time/IrrEnv4/scale1n2.5_2_Res_300Mc5.97e+24_Rc6.37e+06_Lc1e+21k0_0.0458_al-1.68_bt1.87L1e+21_Teq300tau2.txt");//  InitialGuesses/CoredPolytrope/scale1n2.5_2_Res_200Mc5.97e+24_Rc6.37e+06_Lc1e+18k0_151_al-1.00_bt1.00L1e+18_Teq300tau100.txt");//InitialGuesses/CoredPolytrope/scale1n2.5_8.3e+03_Res_400Mc5.97e+24_Rc6.37e+06_Lc1e+20k0_0.00286_al-1.68_bt1.87L1e+20_Ts300Ps5e+03.txt");//"Results/time/IrrEnv3/scale1n2.5_8.3e+03_Res_1000_ftime_7e-310Mc1.19e+25_Rc7.58e+06_Lc3e+19k0_0.00286_al-1.68_bt1.87L3e+19_Teq55tau1t=0.txt"); // ("./"+infolder+"/scale1"+nstring+Resstring+corestring+condstring+Tstring+".txt");
        std::cout << txtstring << std::endl;
	 
    txtfile = &txtstring[0];
    scalefile = &scalestring[0];

    mPr0 = readArray(txtfile);
    scale0 = readArray(scalefile);   
    
    if(J == 100){
	
	planet.y = mPr0.block(0,1,J,I); 
	}else{

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
	}
        planet.scale = scale0.segment(1,I);
	
    ArrayXXd cell_diffs(J-1,I+1);
    cell_diffs.col(0) = planet.m.head(J-1);
    cell_diffs.col(1) = JacobMaker::diff(planet.y.col(0))/planet.y.block(1,0,J-1,1);
    cell_diffs.col(2) = JacobMaker::diff(planet.y.col(1))/planet.y.block(0,1,J-1,1);
    cell_diffs.col(3) = JacobMaker::diff(planet.y.col(2))/planet.y.block(1,2,J-1,1);
    cell_diffs.col(4) = JacobMaker::diff(planet.y.col(3))/planet.y.block(0,3,J-1,1);
    //std::cout << cell_diffs << std::endl;
    
    Allthing.block(0,1,J,I) = planet.y;
    Allthing.col(0) = planet.m;
    
	if ((cell_diffs.block(0,1,J-1,I-1).abs() > 0.25).any() == 1){
        std::cout << "difference between cells larger than 25%"  << std::endl;
	std::cout << Allthing << "\n"<< planet.scale << std::endl;
	//std::cout << cell_diffs << std::endl;

    }
    
	
    planet.L_core = planet.scale(3); // 1400*

    planet.H = 0*planet.scale(3)/planet.M;

    planet.Teq = std::strtod(argv[6] , nullptr ); 
    planet.tau = std::strtod(argv[5] , nullptr );

    if (abs(planet.M - scale0(0)) > planet.M *0.0001){
        std::cout << planet.M << " " << scale0(0) << std::endl;
          throw std::invalid_argument("Wrong mass!" );
    }
    
    double X = 0.05;
    char log_or_lin = '2';
    std::cout << outfolder << std::endl;
    //tStepper.pressure_alter(&planet , 100 , acc , argv[2] , outfolder);
    //tStepper.mu_m_alter(&planet , 2.0 , acc , argv[2] , outfolder);
    //tStepper.n_alter(&planet, 2.5 , acc , argv[2] , outfolder);
    tStepper.timestepper_full(&planet, &EoSIdeal , 1e8 , yrs_max , Lcore_f  ,  acc , ftime , 1 , argv[2] , &log_or_lin , outfolder);
    ///tStepper.tau_alter(&planet , 2.0/3 , acc , argv[2] , outfolder);
    //tStepper.k_alter(&planet , kaim , acc , argv[2] , outfolder);
    //tStepper.L_alter(&planet , 2e16 , acc , argv[2] , outfolder);
    //tStepper.temp_alter(&planet , 800 , acc , argv[2] , outfolder);
    //tStepper.core_incrementer(&planet , 2*M_Earth , X , acc , argv[5] , outfolder );

    cell_diffs.col(0) = planet.m.head(J-1);
    cell_diffs.col(1) = JacobMaker::diff(planet.y.col(0))/planet.y.block(1,0,J-1,1);
    cell_diffs.col(2) = JacobMaker::diff(planet.y.col(1))/planet.y.block(0,1,J-1,1);
    cell_diffs.col(3) = JacobMaker::diff(planet.y.col(2))/planet.y.block(1,2,J-1,1);
    cell_diffs.col(4) = JacobMaker::diff(planet.y.col(3))/planet.y.block(0,3,J-1,1);
    //std::cout << cell_diffs << std::endl;
     Allthing.block(0,1,J,I) = planet.y;
    Allthing.col(0) = planet.m;
	if ((cell_diffs.block(0,1,J-1,I-1).abs() > 0.25).any() == 1){
        std::cout << "difference between cells larger than 25%" << std::endl;
	std::cout << Allthing << "\n"<< planet.scale << std::endl;
    

    }
 
    std::cout << "Completed" << std::endl;
    return 0;
}
