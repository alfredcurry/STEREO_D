#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include "converge.hpp"
#include "Matrix_reader.h"
#include <string>
#include <vector>
#include "Interpolator.hpp"
#include "PhysicalConsts.h"
#include "Jacob_Maker.hpp"
#include "SurfaceTemp.hpp"
#include "SurfacePressure.hpp"
#include "EquationOfState.hpp"

using namespace Eigen;

int main(int arg , char** argv){
    //grid size , mass , EoS type , Tbound type , Pbound type , surface temp OR L, Psurf /opacity, "n" , # 
    int J = std::strtol(argv[1] , nullptr , 0);
    int I = 4;
    double acc = 1e-5;
    std::string folder = "Seager";
    char Tstring[100] , condstring[100] , Resstring[100], nstring[100];
    struct P_P planet(I,J);
    planet.M = std::strtod(argv[2] , nullptr ) * M_Earth;
    planet.M_core = 0;
    planet.R_core = 0;
    planet.L_core = 0;

    planet.Tbound = *argv[4];
    planet.Pbound = *argv[5];

    planet.n = std::strtod(argv[8] , nullptr );         

    double EoSvariables[1];
    EoSvariables[0] = std::strtod(argv[8] , nullptr );
   
    EoS EoSSea(EoSvariables);

    //double kaim = std::strtod(argv[6] , nullptr );
    //planet.conds[1] = std::strtod(argv[7], nullptr);
    //planet.conds[2] = std::strtod(argv[8], nullptr);

    char* num = argv[9];
    char *txtfile , *scalefile;

    ArrayXXd Allthing(J , I + I + 3) , Fls(J,3);

    //m and guesses
    planet.m = ArrayXd::LinSpaced(J,0,1);//::massbins(J);
  
    Converger con(I,J);
    ArrayXXd mPr0;
    ArrayXd scale0; 
    //set up guess
    JacobMaker::prepm( &planet );
  
    //std::cout << J << " " << planet.M <<" " << planet.EoS << " " << planet.Tsurf[0] <<" " << planet.n << std::endl;  
    
    //s - Seagar(2007) (for perovskite) ; I - ideal gas , R = 830
        
    mPr0 = readArray("./InitialGuesses/Seager/results_c_161.txt");
    scale0 = readArray("./InitialGuesses/Seager/scale_c_161.txt");

    ArrayXd mm0 = JacobMaker::diffminus(mPr0.col(0));
    double mtemp = planet.mm(0);
    planet.mm(0) = 0;

    mm0 = mPr0.col(0) - mm0/2 ;
    mm0(0) = 0;
    mm0(mm0.size()-1) = 1;

    planet.y.col(0) = interpolate(planet.mm , mm0 , mPr0.col(1) ); // pressure, using interpolation
    planet.y.block(0,1,J-1,1) = interpolate(planet.m.tail(J-1) , mPr0.col(0) , mPr0.col(2) ); //radius , using interpolation
    planet.mm(0) = mtemp;
 
    planet.y.block(0,3,J-1,1) = planet.m.tail(J-1); 

    double Tc0 , rhoc0;
    int a;

    
        planet.nabtype = 'a';
    planet.y.col(2) = ArrayXd::Constant(J,1.0);
    
        
        if ( planet.Pbound == 'f'){
            planet.Psurf[0] = std::strtod(argv[7] , nullptr ); 
        } else
        if ( planet.Pbound == 'o'){
            planet.opac = std::strtod(argv[7] , nullptr ); //opacity
            planet.conds[0] = 16*Sig_Stefan/(3*planet.opac); 
            planet.conds[1] = -1.0;
            planet.conds[2] = 3.0;
        }
        if ( planet.Tbound == 'f'){
            planet.Tsurf[0] = std::strtod(argv[6] , nullptr );
            planet.scale(2) = planet.Tsurf[0];    

            planet.scale(3) = 1;//(planet.H*planet.dm).sum() * planet.M;
            planet.H = planet.scale(3)/planet.M;
        } else
        if ( planet.Tbound == 'b'){
            planet.scale(3) = std::strtod(argv[6] , nullptr );
            planet.H = planet.scale(3)/planet.M;
        }
        if (planet.Tbound && planet.Pbound == 'f') {
            planet.conds[0] = 0;
            planet.conds[1] = 0;
            planet.conds[2] = 0;
            
            a = snprintf(Tstring, 100 , "Ts%.2g_Ps%.2g" , planet.Tsurf[0] , planet.Psurf[0]);  
             
            planet.scale(0) = planet.Psurf[0]/(planet.y(J-1));
            planet.scale(1) = pow(scale0(1)*pow(scale0(2),4)/planet.scale(0) , 1.0/4);
        }
        if ( planet.Tbound == 'b' && planet.Pbound == 'f'){
            planet.conds[0] = 0;
            planet.conds[1] = 0;
            planet.conds[2] = 0;
            
            planet.scale(0) = pow(planet.scale(3) / (4*M_PI*pow(scale0(2),2 )*Sig_Stefan*pow(Tc0,4)*pow(scale0(1),-0.5)*pow(planet.Psurf[0], 4/(planet.n+1))),1.0/(0.5-4/(planet.n+1))) ;
            planet.scale(1) = pow(scale0(1)*pow(scale0(2),4)/planet.scale(0) , 1.0/4);
            planet.scale(2) = Tc0*scale0(2)/planet.scale(1);    

            a = snprintf(Tstring, 100 , "L%.2g_Ps%.2g" , planet.scale(3), planet.Psurf[0]);  
        }else if (planet.Tbound == 'b' && planet.Pbound == 'o'){
            
            planet.scale(0) = pow(planet.scale(3) / (4*M_PI*pow(scale0(2),2 - 8/(planet.n+1))*Sig_Stefan*pow(Tc0,4)*pow(scale0(1),-0.5-2/(planet.n+1))*pow(2.0/3 * G_Newt*planet.M/planet.opac, 4/(planet.n+1))),1.0/(0.5-4/(planet.n+1))) ;
            planet.scale(1) = pow(scale0(1)*pow(scale0(2),4)/planet.scale(0) , 1.0/4);
            
            planet.scale(2) = Tc0*scale0(2)/planet.scale(1);    

            planet.nabtype = 'a';
            a = snprintf(Tstring, 100 , "L%.2g_opac%.2g" , planet.scale(3), planet.opac);  
        }

        
        if (planet.Tbound == 'b' ){
            planet.Tsurf[0] = T_surface<double>::BB_temp(planet.scale(3)*planet.y(J-2,3) , planet.scale(1)*planet.y(J-2,1));
        }
        if( planet.Pbound == 'o' ){
            planet.Psurf[0] = P_surf<double>::press(planet.M , planet.scale(1) , 0 ,  planet.opac , 0 , 0 , 2.0/3);
        }
        
    
    
    if(planet.nabtype == 'a'){
        a = snprintf(condstring, 100 , "AD" );
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

   //else if (planet.EoS == 'I'){
     //   planet.conds[0] = kaim;
       // planet.rho = Ideal<ArrayXd>::EoSrho(planet.scale(0)*planet.y.col(0) , planet.scale(2)*planet.y.col(2) , planet.boltzrel);
       // planet.k_cond = Conduct<ArrayXd>::k(planet.rho ,planet.scale(2)*planet.y.col(2) , planet.conds[0] , planet.conds[1] ,  planet.conds[2] );
        //while ( planet.k_cond.maxCoeff() > 1e5){
          //  planet.conds[0] = planet.conds[0]/10;
            //planet.rho = Ideal<ArrayXd>::EoSrho(planet.scale(0)*planet.y.col(0) , planet.scale(2)*planet.y.col(2) , planet.boltzrel);
            //planet.k_cond = Conduct<ArrayXd>::k(planet.rho ,planet.scale(2)*planet.y.col(2) , planet.conds[0] , planet.conds[1] ,  planet.conds[2] );
        //}
        //std::cout << "starting k_0 " << planet.conds[0] << std::endl;
    //}
    
    planet.T0 = planet.scale(3)*planet.y.block(0,2,J-1,1);
    planet.P0 = planet.scale(0)*planet.y.block(0,0,J-1,1);

    Map<Matrix<double , Eigen::Dynamic , Eigen::Dynamic , RowMajor> > G(con.Gv, J, I);

    int t = 0;
    double kpres ;
    if (kaim ==0){kpres = 1;}
    else{kpres = kaim/1000;}
    double kstep = 1;
    planet.dt = 1e308;
   
    a = snprintf(Resstring, 100 , "Res_%.1d_" , J);
    
    while (planet.conds[0] <= kaim  + kpres){
        kstep = std::min(kaim - planet.conds[0], planet.conds[0]/100);
        con.nmax = 500;
        con.converge(&planet , &EoSSea, acc );
        
        if (abs(planet.conds[0]-kaim) < kpres || t % 100 ==0 || con.maxinacc > acc){
            //print the stuff
            std::cout << "t: " << t << ", k_cond_0: " << planet.conds[0] << " / " << kaim <<  std::endl;
            if(con.n >  con.nmax){
                std::cout << "not converged to required accuracy in " << con.nmax << std::endl;
            }
            else{
            std::cout << "No. iterations: " << con.n << " of a maximum of: "<< con.nmax << std::endl;
            } 
            
            std::cout << "Max innacuracy: " << con.maxinacc << ", Mean Accuracy: " << con.meaninacc << std::endl; 

            ofstream fout;
            Fls.col(0) = planet.L_conv;
            Fls.col(1) = planet.L_cond;
            Fls.col(2) = planet.ym.col(3);
            fout.open("Fluxes.txt");
                fout << "Fluxes for COND " << planet.conds[0] << "\nL_conv L_cond L_tot\n" << Fls << "\n"; 
            fout.close();

            if (abs(planet.conds[0]-kaim) < kpres || con.maxinacc > 10 * acc ){    
                Allthing.col(0) = planet.m;
                Allthing.block(0,1,J,I) = planet.y;
                Allthing.block(0, I+1,J,1) = planet.rho;
                Allthing.block(0, I+2,J,1) = planet.nabla;
                Allthing.block(0,I+3, J,I) = con.ErrorArray;
                std::cout << Allthing << std::endl;
                if (*num == '0'){
                    std::cout << "no file saved" <<std::endl;
                }else{
                
                    a = snprintf(nstring, 100 , "c%.3g_" , perov_C ) ;
                       
                        
                    fout.open("./Results/"+folder+"/results"+num+nstring+Resstring+condstring+Tstring+".txt");
                    //fout << "m/M $p/P_C$ r/R rho G1 G2 \n";
                        fout << Allthing ;
                    fout.close();
                    fout.open("./Results/"+folder+"/scale"+num+nstring+Resstring+condstring+Tstring+".txt");
                        fout << planet.M <<"\n"<< planet.scale ;
                    fout.close();
                    //fout.open("./TestingRegime/"+folder+"/LTFluxes"+num+nstring+Resstring+condstring+Tstring+".txt");
                      //  fout << Fls;
                    //fout.close();
                }
            }
        }
        if( con.maxinacc > 10*acc){
            std::cout << "too innaccurate\n" << std::endl;
            break;
        }
        else if (con.maxinacc > acc){
            std::cout << "OK to continue\n" << std::endl;
        }
        //"timestep"
        if (kstep < kpres){kstep = 10*kpres;}
        planet.conds[0]= planet.conds[0] + kstep;

        t++; 
    }
    ArrayXXd cell_diffs(J-1,I+1);
    cell_diffs.col(0) = planet.m.head(J-1);
    cell_diffs.block(0,1,J-1,I) = planet.dy/planet.y.block(0,0,J-1,I);

    if ((cell_diffs.block(0,1,J-1,I-1).abs() > 0.25).any() == 1){
        throw std::invalid_argument("too large a difference");
    }
    //std::cout << planet.dm.transpose() << std::endl;
    std::cout << "scales " << planet.scale.transpose() << " Tsurf " << planet.Tsurf[0] << " Psurf " << planet.Psurf[0] << std::endl;

    double Ksurf , Kc;
    Ksurf = pow(planet.Tsurf[0],1+1/planet.n)/pow(planet.Psurf[0],1/planet.n);
    Kc = pow(planet.scale(2),1+1/planet.n)/pow(planet.scale(0),1/planet.n);
    std::cout << "K_surf: " << Ksurf << " K_c: " << Kc << std::endl;
    std::cout << "Kdifference " << abs((Ksurf-Kc)/Ksurf) << std::endl;
    return 0;
}
