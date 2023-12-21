#include <iostream>
#include <Eigen/Dense>
#include "block_tridiag_solve.h"
#include "Jacob_Maker.hpp"
#include "StructStruct.hpp"
#include "converge.hpp"

using namespace Eigen;

void Converger::converge(struct P_P *planet , EoS *eos , double acc){
            std::cout << "converger starting" << std::endl;
            BlockTriDiagSolver<Eigen::Dynamic> Gsolve(J, I); /// Solver of the matrix inversion (written by Richard Booth) which exploits the block tri-diagonal nature. See Bodenheimer.
              
            Map<Matrix<double , Eigen::Dynamic , Eigen::Dynamic , RowMajor> > G(Gv, J, I);
            Map<Matrix<double , Eigen::Dynamic , Eigen::Dynamic , RowMajor> > Gplus(Gvplus, J, I);
            
            double normG , normG0 , f , fplus=0 , slope , alp = 1e-4 , maxinacc0;    
            ArrayXXd divisordy(J,I); 

            /// Calculate the errors for all of the cells, and the squared norm (f)
            f = Jacob.Gvec(planet , eos , Gv );
            
            std::cout << "Psurf calculated: " << planet->Psurf[0] << std::endl;
            std::cout << "Tsurf calculated " << planet->Tsurf[0] << std::endl;  
            std::cout << J << " " << planet->dy.col(0).size() << std::endl;
            Map<Array<double , Eigen::Dynamic , Eigen::Dynamic, RowMajor > > Dy(x, J , I );
            divisordy.block(0,0,J-1,I) = planet->dy;
            divisordy.row(J-1) << planet->y(J-1,0) , 1 , planet->y(J-1,2) , 1;
          
            maxinacc = (G.array()/divisordy).abs().maxCoeff();
            meaninacc = (G.array()/divisordy).abs().mean() ;
            
            maxinacc0 = maxinacc;
            n = 0;
           
            //std::cout << "the G "<< G.array()/divisordy <<"\n"<< std::endl;
            
            Map<Matrix<double , Dynamic , Dynamic , RowMajor> > diffGu(u, J*I, I); //if I wanted to print it
            Map<Matrix<double , Dynamic , Dynamic , RowMajor> > diffGd(d, J*I, I);
            Map<Matrix<double , Dynamic , Dynamic , RowMajor> > diffGl(l, J*I, I);
            
            bool again ;
            again = maxinacc > acc ;
            std::cout << maxinacc << " acc " << acc << "again? " << again << std::endl;

//            std::cout << planet->ym << std::endl;

            while(again == true ) {

                std::cout <<"Iteration: " << n << std::endl;
                std::cout << "scale " << planet->scale.transpose() << " Lcore " << planet->L_core << std::endl;

                Gplus = G; 
                planetplus = *planet;

                if(n !=0){
                    /// If it is not the first step, find and apply a correction using the jacobian calculated on the last step
                    Gsolve.solve(Gvplus, x );
                    if((Dy.abs() > 0.1*planet->y.abs()).any() == 1){
                        std::cout << "Dy reduced as too large relative to y" << std::endl;
                        //std::cout << (Dy/planet->y).transpose() << "\n"<< std::endl; 
                        Dy = Dy / (10 * (Dy / planet->y).eval().abs().maxCoeff());
                    }
                    add_correction(&planetplus , Dy);
                 
                    
                    fplus = Jacob.Gvec(&planetplus , eos , Gvplus );

                    divisordy.block(0,0,J-1,I) = planetplus.dy;
                    divisordy.row(J-1) << planetplus.y(J-1,0) , 1 , planetplus.y(J-1,2) , 1;
                    for(int i=0; i<I ; i++){
                        for(int j=0; j<J-1 ; j++){
                            if(planetplus.dy(j,i) == 0 ){
                            divisordy(j,i) = planetplus.y(j,i)/100;
                            }
                        }
                    }
                    ErrorArray = Gplus.array()/divisordy;
		            maxinacc = ErrorArray.abs().maxCoeff();
                    
                    std::cout << "First try: " << fplus << " " << f << " fplus " << f + alp*slope <<std::endl;
                    std::cout << "maxes " << maxinacc << " " << maxinacc0 << std::endl;

                }
                if(fplus < f + alp * slope && /*maxinacc < maxinacc0 &&  */  n!=0 && n%10 != 0) {
                    std::cout << "accepted restep with same jacobian" << std::endl;
                   
                } 
                else{
                    /// If first step, or old jacobian doesn't work, calculate the jacobian and do the inversion
                    Jacob.GJacob(planet, eos , l , d , u);
                    Gsolve.factor_matrix(l,d,u);
                    
                    Gplus = G;
                    checker(planet , eos , G  , diffGl , acc , int(J-1) , 0 , 0);
			        checker(planet , eos , G  , diffGl , acc , int(J-1) , 0 , 1);
			        checker(planet , eos , G  , diffGl , acc , int(J-1) , 0 , 2);
			        checker(planet , eos , G  , diffGl , acc , int(J-1) , 0 , 3);
                    if(planet->phi(planet->Jcore) < 1.0){
              /*
                    checker(planet , eos , G  , diffGu , acc , int(planet->Jcore+26) , 1 , 0);
			        checker(planet , eos , G  , diffGu , acc , int(planet->Jcore+26) , 1 , 1);
			        checker(planet , eos , G  , diffGu , acc , int(planet->Jcore+26) , 1 , 2);
			        checker(planet , eos , G  , diffGu , acc , int(planet->Jcore+26) , 1 , 3);
                    
                    checker(planet , eos , G  , diffGu , acc , int(planet->Jcore) , 2 , 0);
			        checker(planet , eos , G  , diffGu , acc , int(planet->Jcore) , 2 , 1);
			        checker(planet , eos , G  , diffGu , acc , int(planet->Jcore) , 2 , 2);
			        checker(planet , eos , G  , diffGu , acc , int(planet->Jcore) , 2 , 3);

                    checker(planet , eos , G  , diffGu , acc , int(planet->Jcore) , 1 , 0);
			        checker(planet , eos , G  , diffGu , acc , int(planet->Jcore) , 1 , 1);
			        checker(planet , eos , G  , diffGu , acc , int(planet->Jcore) , 1 , 2);
			        checker(planet , eos , G  , diffGu , acc , int(planet->Jcore) , 1 , 3);
                    
                    
                    checker(planet , eos , G  , diffGu , acc , int(J-2) , 2 , 0);
			        checker(planet , eos , G  , diffGu , acc , int(J-2) , 2 , 1);
			        checker(planet , eos , G  , diffGu , acc , int(J-2) , 2 , 2);
			        checker(planet , eos , G  , diffGu , acc , int(J-2) , 2 , 3);
                    
                    checker(planet , eos , G  , diffGu , acc , int(J-2) , 1 , 0);
			        checker(planet , eos , G  , diffGu , acc , int(J-2) , 1 , 1);
			        checker(planet , eos , G  , diffGu , acc , int(J-2) , 1 , 2);
			        checker(planet , eos , G  , diffGu , acc , int(J-2) , 1 , 3);

                    
                    checker(planet , eos , G  , diffGu , acc , int(planet->Jcore-10) , 1 , 0);
			        checker(planet , eos , G  , diffGu , acc , int(planet->Jcore-10) , 1 , 1);
			        checker(planet , eos , G  , diffGu , acc , int(planet->Jcore-10) , 1 , 2);
			        checker(planet , eos , G  , diffGu , acc , int(planet->Jcore-10) , 1 , 3);
                    
                    checker(planet , eos , G  , diffGu , acc , int(planet->Jcore-10) , 0 , 0);
			        checker(planet , eos , G  , diffGu , acc , int(planet->Jcore-10) , 0 , 1);
			        checker(planet , eos , G  , diffGu , acc , int(planet->Jcore-10) , 0 , 2);
			        checker(planet , eos , G  , diffGu , acc , int(planet->Jcore-10) , 0 , 3);

                    checker(planet , eos , G  , diffGu , acc , int(planet->Jcore+1) , 0 , 0);
			        checker(planet , eos , G  , diffGu , acc , int(planet->Jcore+1) , 0 , 1);
			        checker(planet , eos , G  , diffGu , acc , int(planet->Jcore+1) , 0 , 2);
			        checker(planet , eos , G  , diffGu , acc , int(planet->Jcore+1) , 0 , 3);

                    
*/
                    } 
                    
                    Gsolve.solve(Gvplus, x );
                    //std::cout << Dy.transpose() << std::endl;
                    std::cout << "solve step done" << std::endl;
                    
                    slope = gradient( x , diffGd , diffGu , diffGl , G);

                    if((Dy.abs() > 0.1*planet->y.abs()).any() == 1){
                        std::cout << "Dy reduced as too large relative to y" << std::endl;
                        //std::cout << (Dy/planet->y).transpose() << "\n"<< std::endl; 
                        Dy = Dy / (10 * (Dy / planet->y).eval().abs().maxCoeff());
                    }

                    fplus = f;
                    Gplus = G;
                    std::cout << "starting lsrch" << std::endl;
                    again = lsrch(&planetplus , eos , fplus , Gvplus , slope , Dy , planet);
                    std::cout << "Second try: " << fplus << " " << f << " fplus " << f + alp*slope << std::endl;
                    
                }               
                if( ((f - fplus)/f < 0.001 && n>=nmax) || n > nmax*2){
                    again = false;
                }else if( n >=nmax){ /// If the maximum steps have been exhausted, by the solution seems to be getting better
                    std::cout << "still making progress" << std::endl;
                }
                f = fplus;
                G = Gplus;
                *planet = planetplus;           
                n++;
    
                divisordy.block(0,0,J-1,I) = planet->dy;
                divisordy.row(J-1) << planet->y(J-1,0) , 1 , planet->y(J-1,2) , 1;
                for(int i=0; i<I ; i++){
                    for(int j=0; j<J-1 ; j++){
                        if(planet->dy(j,i) == 0 ){
                            divisordy(j,i) = planet->y(j,i)/100;
                        }
                    }
                }
                
                ErrorArray = G.array()/divisordy;
		        maxinacc = ErrorArray.abs().maxCoeff();
                maxinacc0 = maxinacc;
                if (again == true){ again = (maxinacc > acc);}
                std::cout << maxinacc << " acc " << acc <<std::endl;
            }
            
            meaninacc = ErrorArray.abs().mean() ;

            //if( maxinacc == abs(ErrorArray(0,3))){
              //      std::cout << "Innermost L !!" << std::endl;
                //    maxinacc = ErrorArray.block(1,0,J-1,4).abs().maxCoeff();
            //}
            if ( ( ErrorArray.abs() > acc).any()==1 ){
                if(n<nmax){
                    std::cout<< "Does not want to move" << std::endl; 
                }else{
                    n++;
                }
            }
            if(n > nmax){
                std::cout << "not converged to required accuracy in " << nmax << std::endl;
            }
        
        }

Converger::Converger(int _I , int _J):
            I(_I) , J(_J) , Jacob(_I, _J) , planetplus(_I,_J)
            {
                x = new double[I*J] ;
                l = new double[I*I*J]; 
                d = new double[I*I*J] ;
                u = new double[I*I*J];
                Gv = new double[I*J];
                Gvplus = new double[I*J];
            };

void Converger::add_correction(struct P_P *planet , ArrayXXd Dy ){
            planet->y  =  planet->y - Dy;
        
            planet->scale(0) = planet->scale(0) * planet->y(0,0);
            planet->scale(1) = planet->scale(1) * planet->y(J-1,1);    
            planet->scale(2) = planet->scale(2) * planet->y(0,2);
            planet->scale(3) = planet->scale(3) * planet->y(J-1,3);

            planet->y.col(0) = planet->y.col(0) / planet->y(0,0);
            planet->y.col(1) = planet->y.col(1) / planet->y(J-1,1);
            planet->y.col(2) = planet->y.col(2) / planet->y(0,2);
            planet->y.col(3) = planet->y.col(3) / planet->y(J-1,3);
        }
        
double Converger::gradient(double *x, const MatrixXd &diffGd , const MatrixXd &diffGu , const MatrixXd &diffGl , const MatrixXd &G){
            MatrixXd d_block(I,I) , u_block(I,I) , l_block(I,I);
            RowVectorXd  G_blockd(I) , G_blocku(I) , G_blockl(I) , gradvec(J*I);
            Map<Matrix<double , Eigen::Dynamic , 1 , ColMajor > > Dy(x, J*I , 1 );
            
                d_block = diffGd.block(0,0,I,I);
                l_block = diffGl.block(I,0,I,I);
                G_blockd = G.row(0);
                G_blockl = G.row(1);
                gradvec.segment(0,I) = G_blockd * d_block  + G_blockl * l_block ;
            for( int j=1; j<J-1 ; j++){
                d_block = diffGd.block(j*I,0,I,I);
                u_block = diffGu.block((j-1)*I,0,I,I);
                l_block = diffGl.block((j+1)*I,0,I,I);
                G_blockd = G.row(j);
                G_blocku = G.row(j-1);
                G_blockl = G.row(j+1);
                gradvec.segment(j*I,I) = G_blockd * d_block  + G_blocku * u_block  + G_blockl * l_block ;
            }
                d_block = diffGd.block((J-1)*I,0,I,I);
                u_block = diffGu.block((J-2)*I,0,I,I);
                G_blockd = G.row(J-1);
                G_blocku = G.row(J-2);
                gradvec.segment((J-1)*I,I) = G_blockd * d_block + G_blocku * u_block;
            return -2*gradvec.dot(Dy);
        }
        
bool Converger::lsrch(struct P_P *planetplus , EoS *eos , double &f , double *Gv , const double slope , const ArrayXXd Dy , const struct P_P *planet ){
    /// Algorithm adapted from Numerical Recipes, Press (2007), Third Edition
            double lmda = 1 , lmda2 , f_old , f2 ;
            double TOLX , lmda_tmp , rhs1 , rhs2, a , b , disc , alp = 1e-4;
            TOLX = std::numeric_limits<double>::epsilon();
            int n=0;
            double lmda_min = std::min(1.0e-7 , TOLX*(planet->y/Dy).maxCoeff());
            f_old = f;
            add_correction(planetplus , lmda*Dy);
            f = Jacob.Gvec(planetplus , eos , Gv );
            std::cout << f <<" f? " << f_old  << std::endl;
            while (f > f_old + alp * slope * lmda ){
                planetplus->y = planet->y;
                planetplus->scale = planet->scale;
                if(n == 0){
                    lmda_tmp = -slope/(2.0*(f-f_old-slope));
                } else {
                    rhs1=f-f_old-lmda*slope;
                    rhs2=f2-f_old-lmda2*slope;
                    a=(rhs1/(lmda*lmda)-rhs2/(lmda2*lmda2))/(lmda-lmda2);
                    b=(-lmda2*rhs1/(lmda*lmda)+lmda*rhs2/(lmda2*lmda2))/(lmda-lmda2);
                    if (a == 0.0) lmda_tmp = -slope/(2.0*b);
                    else {
                        disc=b*b-3.0*a*slope;
                    if (disc < 0.0) lmda_tmp=0.5*lmda;
                    else if (b <= 0.0) lmda_tmp=(-b+sqrt(disc))/(3.0*a);
                    else lmda_tmp=-slope/(b+sqrt(disc));
                    }
                    if (lmda_tmp>0.5*lmda)
                        lmda_tmp=0.5*lmda; 
                }
                lmda2 = lmda;
                f2 = f;
                lmda = std::max(lmda_tmp , 0.1*lmda);
                if (lmda <=lmda_min){
                    std::cout<< "believes that mimimum is found (suggested step too small) lmda:" << lmda << " min " << lmda_min << std::endl;
                    add_correction(planetplus , lmda*Dy);
                    f = Jacob.Gvec(planetplus , eos , Gv );
                    return false;
                }
                add_correction(planetplus , lmda*Dy);
                f = Jacob.Gvec(planetplus , eos , Gv );
                n++;
            }
            std::cout << "n " << n << " lamda " << lmda << std::endl;
            return true;
        }
void Converger::checker(struct P_P *planet , EoS *eos , MatrixXd G , MatrixXd diffG , double acc , int row , int Gtype , int col ){
    /* A funcion for checking that derivates are calculated correctly by comparing with brute force finite difference.
        Remember for l should be row-1 and u row+1 (for debugging)
    */
                int row_adjust = -1;
                double smoly = 1e-6*planet->y(row+row_adjust,col);
                double Gv1[I*J] ;
                planet->y(row+row_adjust,col) = planet->y(row+row_adjust,col) + smoly;
                double scale0 = planet->scale(1);
            
            //std::cout << "phi\n" << planet->phi.transpose() << std::endl;
                
                std::cout << "checking " ;
                
                again = Jacob.Gvec(planet , eos , Gv1 );
                Map<Matrix<double , Eigen::Dynamic , Eigen::Dynamic , RowMajor> > G1(Gv1, J, I);
                MatrixXd difG = (G1-G)/smoly;
                //std::cout << difG <<std::endl;
                //std::cout << diffG << std::endl;
                std::cout << "row = " << row << " col = " << col << " type = " << Gtype << " Brute force: " << difG(row,Gtype)<< " calculated: " << diffG((row)*I+Gtype,col) << std::endl;
               
                planet->y(row+row_adjust,col) = planet->y(row+row_adjust,col) - smoly;
        }