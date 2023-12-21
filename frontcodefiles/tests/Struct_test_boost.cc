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

double Gen_Struct( int _J , double n , double Tsurf , double Psurf){
    //grid size , mass , Tbound type , Pbound type , surface temp OR L, Psurf /opacity, "n" , # , boltzrel
    int J = _J;
    int I = 4;
    double acc = 1e-5;
    std::string folder ;
    folder = "Tests/RunFiles";

    char Tstring[100] , condstring[100] , Resstring[100], nstring[100];
    struct P_P planet(I,J);
    planet.M = M_Earth;
    planet.M_core = 0;
    planet.R_core = 0;
    
    planet.Tbound = 'f';
    planet.Pbound = 'f';
    planet.n = n;      
    planet.mu_m = 1.0;

    planet.BL = '0';
    planet.MassLossType = '0';

    ArrayXXd Allthing(J , I + I + 3) , Fls(J,3);

    double EoSvariables[2];
    EoSvariables[0] = planet.mu_m ;
    EoSvariables[1] = planet.n ;
    EoS EoSIdeal(EoSvariables);
    folder.append(EoSIdeal.EoStype);

    //m and guesses
    planet.m = JacobMaker::massbins(J); // ArrayXd::LinSpaced(J,0,1);//
  
    Converger con(I,J);
    ArrayXXd mPr0;
    ArrayXd scale0; 
  
    //set up guess
    JacobMaker::prepm( &planet );
  
    std::string txtstring = ("./InitialGuesses/ColdPolytrope/results1000_n_"+std::to_string(static_cast<int>(planet.n*10+0.5))+".txt");
    std::string scalestring = ("./InitialGuesses/ColdPolytrope/scale1000_n_"+std::to_string(static_cast<int>(planet.n*10+0.5))+".txt");
    std::cout << txtstring << std::endl;
    char *txtfile = &txtstring[0];
    char *scalefile = &scalestring[0];

    mPr0 = readArray(txtfile);
    scale0 = readArray(scalefile);
    
    ArrayXd mm0 = JacobMaker::diffminus(mPr0.col(0));
    double mtemp = planet.mm(0);
    planet.mm(0) = 0;

    mm0 = mPr0.col(0) - mm0/2 ;
    mm0(0) = 0;
    mm0(mm0.size()-1) = 1;

    planet.y.col(0) = interpolate(0.98*planet.mm , mm0 , mPr0.col(1) ); // pressure, using interpolation
    planet.y.block(0,1,J-1,1) = interpolate(0.98*planet.m.tail(J-1) , mPr0.col(0) , mPr0.col(2) ); //radius , using interpolation
    planet.mm(0) = mtemp;
 
    planet.y.block(0,3,J-1,1) = planet.m.tail(J-1); 

    double Tc0 , rhoc0;
    int a;
    
    planet.nabtype = 'a';

    planet.y.col(2) = pow(planet.y.col(0),1.0/(planet.n+1));
        
    planet.Psurf[0] = Psurf;     
    planet.Tsurf[0] = Tsurf;
    planet.scale(3) = 1;//(planet.H*planet.dm).sum() * planet.M;
    planet.H = planet.scale(3)/planet.M;
        
    rhoc0 = mPr0(0,3)*scale0(0)/pow(scale0(2),3);
    Tc0 = scale0(1)/(rhoc0*boltz_R/planet.mu_m);
    
    planet.conds[0] = 0;
    planet.conds[1] = 0;
    planet.conds[2] = 0;
             
    planet.scale(0) = planet.Psurf[0]/(planet.y(J-1));
    planet.scale(1) = pow(scale0(1)*pow(scale0(2),4)/planet.scale(0) , 1.0/4);
    planet.scale(2) = Tc0*scale0(2)/planet.scale(1);    

    a = snprintf(nstring, 100 , "n%.1f_%.3g_" , planet.n , boltz_R/planet.mu_m) ;
    a = snprintf(Tstring, 100 , "Ts%.2g_Ps%.2g" , planet.Tsurf[0] , planet.Psurf[0]); 
    a = snprintf(condstring, 100 , "AD" );
    a = snprintf(Resstring, 100 , "Res_%.1d_" , J);
    
    planet.y(J-1,1) = planet.y(J-2,1);
    planet.y(J-1,3) = planet.y(J-2,3);
    planet.y(J-1,0) = planet.Psurf[0]/planet.scale(0);
    planet.y(J-1,2) = planet.Tsurf[0]/planet.scale(2);
    
    double kaim = planet.conds[0];
    std::cout << "\n" << planet.scale << std::endl;
    std::cout << planet.Tsurf[0] << " T P " << planet.Psurf[0] << std::endl;
    planet.nu = ArrayXd::Constant(J,5.0e11);
    
    planet.T0 = planet.scale(3)*planet.y.block(0,2,J-1,1);
    planet.P0 = planet.scale(0)*planet.y.block(0,0,J-1,1);

    Map<Matrix<double , Eigen::Dynamic , Eigen::Dynamic , RowMajor> > G(con.Gv, J, I);

    planet.dt = 1e308;
   
      
    con.nmax = 500;
    con.converge(&planet , &EoSIdeal, acc );
                
    //print the stuff
    if(con.n >  con.nmax){
          std::cout << "not converged to required accuracy in " << con.nmax << std::endl;
    }
    else{
          std::cout << "No. iterations: " << con.n << " of a maximum of: "<< con.nmax << std::endl;
    
            
          std::cout << "Max innacuracy: " << con.maxinacc << ", Mean Accuracy: " << con.meaninacc << std::endl; 

          ofstream fout;
              
          Allthing.col(0) = planet.m;
          Allthing.block(0,1,J,I) = planet.y;
          Allthing.block(0, I+1,J,1) = planet.rho;
          Allthing.block(0, I+2,J,1) = planet.nabla;
          Allthing.block(0,I+3, J,I) = con.ErrorArray;
                        
          fout.open("./"+folder+"/resultsSTRUCT_TEST"+nstring+Resstring+condstring+Tstring+".txt");
                //fout << "m/M $p/P_C$ r/R rho G1 G2 \n";
                fout << Allthing ;
                fout.close();
                fout.open("./"+folder+"/scaleSTRUCT_TEST"+nstring+Resstring+condstring+Tstring+".txt");
                fout << planet.M <<"\n"<< planet.scale ;
                fout.close();
               
    }
    if( con.maxinacc > acc){
            std::cout << "too innaccurate\n" << std::endl;
    }
    
    std::cout << "scales " << planet.scale.transpose() << " Tsurf " << planet.Tsurf[0] << " Psurf " << planet.Psurf[0] << std::endl;

    ArrayXd rho_th ;
    double K0 = std::pow(planet.Tsurf[0]*boltz_R/planet.mu_m ,1+1.0/planet.n)/std::pow(planet.Psurf[0] , 1.0/planet.n);
    double radius_const = sqrt(2*M_PI*G_Newt/K0);
    std::cout << "radius const " <<  radius_const << std::endl;
    rho_th = planet.rho(0) * sin(radius_const*planet.ym.col(1))/(radius_const*planet.ym.col(1));
    std::cout << "Static_Poly_Test completed\n" << std::endl;
    
    return sqrt(pow((rho_th - planet.rho)/rho_th,2.0).sum()/J);
}

BOOST_AUTO_TEST_CASE( Static_Poly_Test )
{
  boost::unit_test::unit_test_log.set_stream( std::cerr );
  boost::unit_test::unit_test_log.set_threshold_level(boost::unit_test::log_level::log_test_units); 
//  namespace tt = boost::test_tools;

  BOOST_TEST(Gen_Struct(500 , 1 , 10 , 1e5 ) < 1.7e-2);
  // For reference:
  //[  2000 , 1500 ,  1000 , 500 , 200 ,  100 , 10  ]
  //[3.86409282e-03 4.91592581e-03 7.17564442e-03 1.60787266e-02 7.02595918e-02 3.50311173e-01 4.37773709e+00]
}
