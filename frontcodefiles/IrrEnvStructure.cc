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
    //grid size , mass , surface temp ,  surface P , # , tmax , Mcore , Rcore , L 
    int J = std::strtol(argv[1] , nullptr , 0);
    int I = 4;
    int a;
    double acc = 1e-5 ;
    double Taim = 375 , Paim = 35 , tmax = std::strtol(argv[6] , nullptr , 0) , Laim = 1e22;
    std::string infolder , outfolder ; 
    char *txtfile , *scalefile;

    outfolder =  "Results/time/IrrEnv";
    infolder = "InitialGuesses/CoredPolytrope";
    ArrayXXd mPr0;
    ArrayXd scale0; 
    struct P_P planet(I,J);
    stepper tStepper(I,J);

    //planet properties
    planet.M = std::strtod(argv[2] , nullptr ) * M_Earth;
    
    planet.n = 1.5;
    planet.mu_m = 1;
    planet.M_core = std::strtod(argv[7] , nullptr )*M_Earth;
    planet.R_core = std::strtod(argv[8] , nullptr )*R_Earth;
    planet.scale(3) = std::strtod(argv[9] , nullptr ); 
	planet.L_core = planet.scale(3);
    double alpha_Owen , beta_Owen , kap0_Owen;
    alpha_Owen = 0.68;
    beta_Owen = 0.45;
    kap0_Owen = 2.29e-7; //in SI : cgs - 0.1Pa 4.79e-8*10^alpha?
    planet.Tbound = 'f';
    planet.Pbound = 'f';
    planet.nabtype = 'c';

    planet.nu = ArrayXd::Constant(J,1.0e-18);

    double EoSvariables[2];
    EoSvariables[0] = planet.mu_m ;
    EoSvariables[1] = planet.n ;
    EoS EoSIdeal(EoSvariables);

    planet.Psurf[0] = std::strtod(argv[4] , nullptr );
    planet.Tsurf[0] = std::strtod(argv[3] , nullptr );

    double kaim = 16*Sig_Stefan/(3 * kap0_Owen) * pow(boltz_R/planet.mu_m, -alpha_Owen);
    planet.conds[0] = kaim ;    
    planet.conds[1] = -1-alpha_Owen;        
    planet.conds[2] = 3-alpha_Owen-beta_Owen;
    //read in guesses
    std::string txtstring , scalestring;
    if (planet.M_core ==0 && planet.Psurf[0] == 1e7){
        txtstring = ("./InitialGuesses/ColdPolytrope/results1000_n_"+std::to_string(static_cast<int>(planet.n*10+0.5))+".txt");
        scalestring= ("./InitialGuesses/ColdPolytrope/scale1000_n_"+std::to_string(static_cast<int>(planet.n*10+0.5))+".txt");
    
    }else 
    if(planet.M_core ==0 && planet.Psurf[0] == 100){
        char Tstring[100] , condstring[100] , nstring[100] , Resstring[100] , corestring[100];
        a = snprintf(nstring, 100 , "n%.1f_%.3g_" , planet.n , boltz_R/planet.mu_m) ;
        a = snprintf(Tstring, 100 , "L%.2g_Ts%.3gPs%.2g" , planet.scale(3) , planet.Tsurf[0] , planet.Psurf[0] );  
        a = snprintf(corestring, 100 , "_" );  
	if(planet.nabtype == 'a'){
                a = snprintf(condstring, 100 , "AD" );
        }else{
                a = snprintf(condstring, 100 , "k0_%.3g_al%.2f_bt%.2f" , planet.conds[0], planet.conds[1] , planet.conds[2]);
        }
        a = snprintf(Resstring, 100 , "Res_%.1d" , J);
        txtstring = ("./"+infolder+"/results2"+nstring+Resstring+corestring+condstring+Tstring+".txt");
        scalestring = ("./"+infolder+"/scale2"+nstring+Resstring+corestring+condstring+Tstring+".txt");
        std::cout << txtstring << std::endl;
    }else if(*argv[5] == '5'){	
	txtstring = "InitialGuesses/CoredPolytrope/results1n1.5_8.3e+03_Res_500Mc5.97e+24_Rc6.37e+06k0_0.00286_al-1.68_bt1.87Ts1.5e+02_Ps1.5e+04_L1e+16.txt";
	scalestring = "InitialGuesses/CoredPolytrope/scale1n1.5_8.3e+03_Res_500Mc5.97e+24_Rc6.37e+06k0_0.00286_al-1.68_bt1.87Ts1.5e+02_Ps1.5e+04_L1e+16.txt";
	}else{
        char Tstring[100] , condstring[100] , nstring[100] , Resstring[100] , corestring[100];
        a = snprintf(nstring, 100 , "n%.1f_%.3g_" , planet.n , boltz_R/planet.mu_m) ;
        a = snprintf(Tstring, 100 , "L%.2g_Ts%.3gPs%.2g"  , planet.scale(3) , planet.Tsurf[0] , planet.Psurf[0]); 
        a = snprintf(corestring, 100 , "Mc%.3g_Rc%.3g_Lc%.2g" , planet.M_core , planet.R_core , planet.L_core);  
        if(planet.nabtype == 'a'){
		a = snprintf(condstring, 100 , "AD" );
	}else{
		a = snprintf(condstring, 100 , "k0_%.3g_al%.2f_bt%.2f" , planet.conds[0], planet.conds[1] , planet.conds[2]);
        }
	a = snprintf(Resstring, 100 , "Res_%.1d" , J);
        txtstring = ("./"+infolder+"/results3"+nstring+Resstring+corestring+condstring+Tstring+".txt");
        scalestring = ("./"+infolder+"/scale3"+nstring+Resstring+corestring+condstring+Tstring+".txt");
        std::cout << txtstring << std::endl;
    }
    txtfile = &txtstring[0];
    scalefile = &scalestring[0];

    mPr0 = readArray(txtfile);
    std::cout << mPr0.col(0).size() << " " << mPr0(mPr0.col(0).size()-1,0) << std::endl;
	scale0 = readArray(scalefile);   
	
	planet.nabtype = 'c';
	planet.L_core = planet.scale(3);
	
    if (planet.M_core ==0 && planet.Psurf[0] == 1e7){
       

        planet.m = ArrayXd::LinSpaced(J,0,1);//JacobMaker::massbins(J);
        JacobMaker::prepm( &planet );
        ArrayXd mm0 = JacobMaker::diffminus(mPr0.col(0));

        double mtemp = planet.mm(0);
        planet.mm(0) = 0;

        mm0 = mPr0.col(0) - mm0/2 ;
        mm0(0) = 0;
        mm0(mm0.size()-1) = 1;
        
        planet.y.col(0) = interpolate(0.975*planet.mm , mm0 , mPr0.col(1) ); // pressure, using interpolation
        planet.y.block(0,1,J-1,1) = interpolate(0.975*planet.m.tail(J-1) , mPr0.col(0) , mPr0.col(2) ); //radius , using interpolation
        //planet.y.block(0,3,J-1,1) = planet.m.tail(J-1); 
        planet.mm(0) = mtemp;
    
        double Tc0 , rhoc0;     
        rhoc0 = mPr0(0,3)*scale0(0)/pow(scale0(2),3);
        Tc0 = scale0(1)/(rhoc0*boltz_R/planet.mu_m);
        planet.scale(0) = planet.Psurf[0]/(planet.y(J-1,0));

        planet.scale(1) = pow(scale0(1)*pow(scale0(2),4)/planet.scale(0) , 1.0/4);
        planet.scale(2) = Tc0*scale0(2)/planet.scale(1);    
        

        planet.y.col(2) = pow(planet.y.col(0),1.0/(planet.n+1));
        planet.y(J-1,1) = planet.y(J-2,1);
        //planet.y(J-1,3) = planet.y(J-2,3);
        planet.y(J-1,0) = planet.Psurf[0]/planet.scale(0);
        planet.y(J-1,2) = planet.Tsurf[0]/planet.scale(2);     
        planet.y.col(3) = ArrayXd::Constant(J,1);
    }else{
//	int J_bound = 50;
  //  	double m_bound = 0.001;
//	planet.m.head(J_bound) = ArrayXd::LinSpaced(J_bound,planet.M_core/planet.M,planet.M_core/planet.M+m_bound);
  //  	planet.m.tail(J-J_bound) = ArrayXd::LinSpaced(J-J_bound,planet.M_core/planet.M+m_bound + m_bound/J_bound, 1);
           
    //	JacobMaker::prepm( &planet );

    //	ArrayXd mm0 = JacobMaker::diffminus(mPr0.col(0));

      //  double mtemp = planet.mm(0);
       
      //  planet.mm(0) = planet.M_core/planet.M;
       // mm0 = mPr0.col(0) - mm0/2 ;
      //  mm0(0) = planet.mm(0);
       // mm0(mm0.size()-1) = 1;
        
       // planet.y.col(0) = interpolate(planet.mm , mm0 , mPr0.col(1) ); // pressure, using interpolation
       // planet.y.block(0,1,J-1,1) = interpolate(planet.m.tail(J-1) , mPr0.col(0) , mPr0.col(2) ); //radius , using interpolation
       // planet.y(J-1,1) = planet.y(J-2,1);
        //planet.y.col(2) = interpolate(planet.mm , mm0 , mPr0.col(3) ); // temp, using interpolation
        //planet.y.block(0,3,J-1,1) = interpolate(planet.m.tail(J-1) , mPr0.col(0) , mPr0.col(4) ); //L, using interpolation
        //planet.y(J-1,3) = planet.y(J-2,3);
       
        //planet.mm(0) = mtemp;
        planet.scale = scale0.segment(1,I);
        planet.y = mPr0.block(0,1,J,I);

std::cout << "here " << std::endl;
        planet.m = mPr0.col(0);
std::cout << mPr0.col(0).size() << std::endl;
        JacobMaker::prepm( &planet );
    }
     std::cout << planet.scale << "\n" << planet.y << std::endl;
    
    
    planet.H = 0*planet.scale(3)/planet.M;
    
    if (abs(planet.M - scale0(0)) > planet.M *0.0001){
        std::cout << planet.M << " " << scale0(0) << std::endl;
          throw std::invalid_argument("Wrong mass!" );
    }

    double X = 0.1;
   
    //tStepper.temp_alter(&planet , Taim , acc , argv[5] , outfolder);
    
    //tStepper.pressure_alter(&planet , Paim , acc , argv[5] , outfolder);
    //tStepper.core_incrementer(&planet , M_Earth , X , acc , argv[5] , outfolder );
    tStepper.L_alter(&planet , &EoSIdeal , Laim , acc , argv[5] , outfolder);
    
	std::cout << "Completed" << std::endl;
    return 0;
}
