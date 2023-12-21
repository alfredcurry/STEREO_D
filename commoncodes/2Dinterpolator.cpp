#include "Matrix_reader.h"
#include <stdexcept>
#include <string>
#include <fstream>
#include "2Dinterpolator.hpp"

template<class Real>
    void Interpolators::BicubicInterpolator<Real>::save_coeff(std::string &filebase, std::string &F_name){
        std::ofstream fout;
        f_name = F_name;
        fout.open(filebase+"/coeff/"+f_name+std::to_string(X.size())+"x"+std::to_string(Y.size())+".txt");
            for(int i=0; i<_a.size(); i++){
                fout << _a[i] << std::endl;
            }
        fout.close();
    
    }

template<class Real>
    void Interpolators::BicubicInterpolator<Real>::read_setup(std::string &filebase , std::string &x_name , std::string &y_name , std::string &F_name , int *grid_size){

        X_name = x_name;
        Y_name = y_name;
        f_name = F_name;
        read_XY(filebase , grid_size);
        Eigen::ArrayXXd f  , f_dX , f_dY , f_dXdY ;
        std::string filename;

        filename = filebase+"/data/"+f_name+std::to_string(grid_size[0])+"x"+std::to_string(grid_size[1])+".txt";
        std::cout << "Setup read from:\n" << filename << std::endl;
        f = readArray(&filename[0]);

        filename = filebase+"/data/"+f_name+"_d"+X_name+std::to_string(grid_size[0])+"x"+std::to_string(grid_size[1])+".txt";
        std::cout << "derivs from:\n" << filename << std::endl;
        
        f_dX = readArray(&filename[0]);
        
        filename = filebase+"/data/"+f_name+"_d"+Y_name+std::to_string(grid_size[0])+"x"+std::to_string(grid_size[1])+".txt";
        f_dY = readArray(&filename[0]);

        filename = filebase+"/data/"+F_name+"_d"+X_name+"d"+Y_name+std::to_string(grid_size[0])+"x"+std::to_string(grid_size[1])+".txt";
        f_dXdY = readArray(&filename[0]);
        setupInterpolator( f , f_dX , f_dY , f_dXdY );

    }
  
template<class Real>
    void Interpolators::BicubicInterpolator<Real>::setupInterpolator(std::vector<Real> x , std::vector<Real> y , Eigen::ArrayXXd f, Eigen::ArrayXXd f_dx , Eigen::ArrayXXd f_dy , Eigen::ArrayXXd f_dxdy)
    {
        X = x;
        Y = y;
        
        setupInterpolator(f , f_dx , f_dy , f_dxdy);
    }
    
template<class Real>
    void Interpolators::BicubicInterpolator<Real>::setupInterpolator(Eigen::ArrayXXd f, Eigen::ArrayXXd f_dx , Eigen::ArrayXXd f_dy , Eigen::ArrayXXd f_dxdy)
    {
        std::cout << "ready to set up" << std::endl;

        X_range[0] = *std::min_element(X.begin(),X.end());
        X_range[1] = *std::max_element(X.begin(),X.end());
        Y_range[0] = *std::min_element(Y.begin(),Y.end());
        Y_range[1] = *std::max_element(Y.begin(),Y.end());

        _a.resize(16*(X.size()-1)*(Y.size()-1));
        Real *a_data = &_a.front();

        //_a[16*X.size()*Y.size()];
        
    
        Matrix44 Left, Right;
        Real xL ,yL ; 
        Matrix44 F;

        Left <<  1,  0,  0,  0,
                 0,  0,  1,  0,
                -3,  3, -2, -1,
                 2, -2,  1,  1;

        Right<<  1,  0, -3,  2,
                 0,  0,  3, -2,
                 0,  1, -2,  1,
                 0,  0, -1,  1;

        for(int i = 0; i < X.size() - 1; i++){
            for( int j = 0; j < Y.size() - 1; j++){

                Eigen::Map<Matrix44>a(a_data + 16*(Y.size()-1)*i + 16*j, 4, 4) ;
                a = Left; 

                xL = X[i+1] - X[i];
                yL = Y[j+1] - Y[j];
                F.block(0,0,2,2) = f.block(i,j,2,2);
                F.block(2,0,2,2) = xL*f_dx.block(i,j,2,2);
                F.block(0,2,2,2) = yL*f_dy.block(i,j,2,2);
                F.block(2,2,2,2) = xL*yL*f_dxdy.block(i,j,2,2);

                a = Left * F * Right;
    
            }
        }
    
    }
