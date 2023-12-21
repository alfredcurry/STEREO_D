#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include "converge.h"
#include "Matrix_reader.h"

using namespace Eigen;

int main(){
    int J = 100, I = 3;
    double acc = 1e-5;

    struct P_P planet;

    ArrayXXd Allthing(J , I + I + 2) ;
    char *txtfile , *scalefile;
    //m and guesses
    planet.m = ArrayXd::LinSpaced(J, 1.0/J , 1);

    planet.EoS = 'I'; 
    planet.boltzrel = 8300;

    ArrayXd Tvec(5);
    Tvec << 10 , 50, 100 , 500 , 1000;

    planet.y = ArrayXXd::Zero(J,I);
    planet.dy = ArrayXXd::Zero(J,I);
    planet.scale = ArrayXd::Zero(I);

    Converger con(I,J);            

    //"timesteping"
    Map<Matrix<double , Eigen::Dynamic , Eigen::Dynamic , RowMajor> > G(con.Gv, J, I);
    ofstream fout;

    double nstep = 0.1;

    for (int i = 0; i < Tvec.size(); i++){ 
        planet.n = 1; 
        planet.Tsurf = Tvec(i);
        
        for (int t = 0; t<10 ; t++){
        
            //set up guess
            con.prepm( &planet );
            std::string txtstring = ("../Cold_Struct/Results/Polytrope/results_n_"+std::to_string(static_cast<int>(planet.n*10+0.5))+".txt");
            std::string scalestring = ("../Cold_Struct/Results/Polytrope/scale_n_"+std::to_string(static_cast<int>(planet.n*10+0.5))+".txt");

            txtfile = &txtstring[0];
            scalefile = &scalestring[0];

            ArrayXXd mPr0 = readArray(txtfile);
            ArrayXd scale0 = readArray(scalefile);

            planet.y.col(0) = mPr0.col(1); //guess pressure
            planet.y.col(1) = mPr0.col(2); //guess radius
            planet.y.col(2) = pow(planet.y.col(0),1/(planet.n+1));

            planet.Psurf = scale0(1)*planet.y(J-1,0)*pow(planet.Tsurf * planet.boltzrel * mPr0(J-1,3)*M_Earth/ (pow(scale0(2),3)*scale0(1)*planet.y(J-1,0)) , 4); 

            planet.scale(2) = planet.Tsurf/planet.y(J-1,2); 
            planet.scale(0) = planet.Psurf/planet.y(J-1,0);
    
            planet.scale(1) = pow(M_Earth*mPr0(J-1,3)/( planet.Psurf / (planet.Tsurf * planet.boltzrel) ) , 1.0/3);  

            planet.T0 = planet.y.col(2);
            planet.P0 = planet.y.col(0);

            con.converge(&planet , acc );

            //print the stuff
            std::cout << "t: " << t << ", polytrope: " << 10*planet.n << " " << static_cast<int>(10*planet.n+0.5) << std::endl;
            if(con.n >= con.nmax){
                std::cout << "not converged to required accuracy in " << con.n << std::endl;
            }
            else{
            std::cout << "No. iterations: " << con.n << std::endl;
            } 
            std::cout << "Max innacuracy: " << (G.array()/planet.dy).abs().maxCoeff() << ", Mean Accuracy: " <<(G.array()/planet.dy).abs().mean() << "\n" << std::endl; 

            Allthing.col(0) = planet.m;
            Allthing.block(0,1,J,I) = planet.y;
            Allthing.block(0, I+1,J,1) = planet.rho;
            Allthing.block(0,I+2, J,I) = (G.array()/planet.dy) ;
   
            fout.open("./results/Polytrope/DirT_"+std::to_string(static_cast<int>(planet.Tsurf))+"_results_n_"+std::to_string(static_cast<int>(planet.n*10+0.5))+".txt");
                //fout << "m/M $p/P_C$ r/R rho G1 G2 \n";
                fout << Allthing ;
            fout.close();
            fout.open("./results/Polytrope/DirT_"+std::to_string(static_cast<int>(planet.Tsurf))+"_scale_n_"+std::to_string(static_cast<int>(planet.n*10+0.5))+".txt");
                fout << M_Earth <<"\n"<< planet.scale ;
            fout.close();
            
        
            //"timestep"
            planet.n = planet.n + nstep;
        
        }
    }
    return 0;
}