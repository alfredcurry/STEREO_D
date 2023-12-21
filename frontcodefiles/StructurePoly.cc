#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include "converge.h"
#include "Matrix_reader.h"
#include "PrepPlanet.h"

using namespace Eigen;

int main(){
    int J = 100, I = 4;
    double acc = 1e-5;

    struct P_P planet(I,J);

    ArrayXXd Allthing(J , I + I + 2) ;

    //m and guesses
    planet.m = ArrayXd::LinSpaced(J, 1.0/J , 1);

    planet.EoS = 'I';
    planet.n = 1.0;
    planet.M = M_Earth;

    Converger con(I,J);
    
    //set up guess
    Prep_Planet::prepm( &planet );

    char *txtfile , *scalefile;
    
    std::string txtstring = ("../Cold_Struct/Results/Polytrope/results_n_"+std::to_string(static_cast<int>(planet.n*10+0.5))+".txt");
    std::string scalestring = ("../Cold_Struct/Results/Polytrope/scale_n_"+std::to_string(static_cast<int>(planet.n*10+0.5))+".txt");

    txtfile = &txtstring[0];
    scalefile = &scalestring[0];

    ArrayXXd mPr0 = readArray(txtfile);
    ArrayXd scale0 = readArray(scalefile);
    
    planet.y.col(0) = mPr0.col(1); //guess pressure
    planet.y.col(1) = mPr0.col(2); //guess radius
    planet.y.col(2) = pow(planet.y.col(0),1/(planet.n+1));
    planet.y.col(3) = planet.m; 
    planet.boltzrel = 8300;
    planet.Tsurf = 100;
    planet.Psurf = scale0(1)*planet.y(J-1,0)*pow(planet.Tsurf * planet.boltzrel * mPr0(J-1,3)*planet.M/ (pow(scale0(2),3)*scale0(1)*planet.y(J-1,0)) , 4); 

    planet.scale(2) = planet.Tsurf/planet.y(J-1,2); 
    planet.scale(0) = planet.Psurf/planet.y(J-1,0);
    planet.scale(1) = pow(planet.M*mPr0(J-1,3)/( planet.Psurf / (planet.Tsurf * planet.boltzrel) ) , 1.0/3);              
    planet.scale(3) = 4*M_PI*pow(planet.scale(1),2) * Sig_Stefan * pow(planet.Tsurf,4) / planet.y(J-1,3);

//    planet.scale(1) = sqrt(M_PI /(2*G_Newt)*pow(planet.boltzrel*planet.scale(2)/planet.scale(0),1+1/planet.n)*planet.scale(0));
    planet.H = planet.scale(3)/(J*planet.M*planet.dm);
    planet.nu = ArrayXd::Constant(J,5.0e11);
    planet.T0 = planet.y.col(2);
    planet.P0 = planet.y.col(0);
    planet.kap = ArrayXd::Constant(J,0.01*3.0/(3000*1000)); 
    //"timesteping"
    Map<Matrix<double , Eigen::Dynamic , Eigen::Dynamic , RowMajor> > G(con.Gv, J, I);
    ofstream fout;

    double nstep = 0.01;
    planet.dt = 3600.0*24*365*10000;
    for (int t = 0; t<1 ; t++){
        
        con.converge(&planet , acc );

        if (t % 10 == 0){
            //print the stuff
            std::cout << "t: " << t << ", polytrope: " << 10*planet.n << " " << static_cast<int>(10*planet.n+0.5) << std::endl;
            if(con.n >= con.nmax){
                std::cout << "not converged to required accuracy in " << con.n << std::endl;
            }
            else{
            std::cout << "No. iterations: " << con.n << " of a maximum of: "<< con.nmax << std::endl;
            } 
            std::cout << "Max innacuracy: " << (G.array()/planet.dy).abs().maxCoeff() << ", Mean Accuracy: " <<(G.array()/planet.dy).abs().mean() << "\n" << std::endl; 

            Allthing.col(0) = planet.m;
            Allthing.block(0,1,J,I) = planet.y;
            Allthing.block(0, I+1,J,1) = planet.rho;
            Allthing.block(0,I+2, J,I) = (G.array()/planet.dy) ;
   
            fout.open("./results/Polytrope/results3_n_"+std::to_string(static_cast<int>(planet.n*10+0.5))+".txt");
                //fout << "m/M $p/P_C$ r/R rho G1 G2 \n";
                fout << Allthing ;
            fout.close();
            fout.open("./results/Polytrope/scale3_n_"+std::to_string(static_cast<int>(planet.n*10+0.5))+".txt");
                fout << planet.M <<"\n"<< planet.scale ;
            fout.close();
            
        }
        //"timestep"
        planet.n = planet.n + nstep;
        
    }

    return 0;
}