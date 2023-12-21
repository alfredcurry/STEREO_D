#ifndef Interpolators_BicubicInterpolator_hpp
#define Interpolators_BicubicInterpolator_hpp

#include <iostream>
#include <Eigen/Dense>
#include <stdexcept>
#include <string>
#include <fstream>
#include <vector>
#include "Vector_reader.h"
#include "File_error.h"

namespace Interpolators{

    template<class Real> /// Interpolation on a 2D grid, where values and derivatives are known
    class BicubicInterpolator{
        public: 

            using Matrix44 = Eigen::Matrix<Real,4,4 >;
            using Matrix22 = Eigen::Matrix<Real,2,2 >;
            using Matrix44Array = Eigen::Array< Matrix44, Eigen::Dynamic, Eigen::Dynamic >;
            using ColVector4 = Eigen::Matrix<Real,4,1 >;
            using RowVector4 = Eigen::Matrix<Real,1,4 >;

            std::vector<Real> X , Y;
            /// @brief Setup interpolator coefficients from inputted arrays
            void setupInterpolator(std::vector<Real> x ,std::vector<Real> y , Eigen::ArrayXXd f , Eigen::ArrayXXd f_dx , Eigen::ArrayXXd f_dy , Eigen::ArrayXXd f_dxdy);
            int* return_index(Real x , Real y ); /// Return the grid indices closest to x,y
            Real return_value(Real x , Real y , int i , int j); /// Return the interpolated value at x,y
            Real* return_grads(Real x , Real y , int i , int j); /// Return the interpolated gradient (returns array of size 2) at x,y
            void read_setup(std::string &filebase,std::string &x_name , std::string &y_name , std::string &F_name , int *grid_size); /// Setup interpolator from read files
            void save_coeff(std::string &filebase, std::string &F_name); /// Save interpolator coefficients

            void read_coeff(std::string &coeffile);  /// Read interpolator coefficients from a file
            void setup_from_coeff(std::string &filebase,std::string &x_name , std::string &y_name , std::string &F_name , int *grid_size); /// Setup the interpolator by reading the coefficients from a file
       
        private:
            void read_XY(std::string &filebase, int *grid_size); /// Read the grid x,y from a file
            void setupInterpolator(Eigen::ArrayXXd f , Eigen::ArrayXXd f_dx , Eigen::ArrayXXd f_dy , Eigen::ArrayXXd f_dxdy); ///setup the interpolator coefficients from arrays (grid already known)

            int get_x_index_to_left_of(Real x) const /// Find the closest x-index of grid to x
            {
                //assume regular grid for now
                int index = (x-X_range[0])/(X_range[1]-X_range[0])*(X.size()-1);
                if(index == (X.size()-1)){
                    return index - 1;
                }
                else{
                    return index;
                } 
            }
            int get_y_index_below(Real y) const /// Find the closest y-index of grid to y
            {
                //assume regular grid for now
                int index = (y-Y_range[0])/(Y_range[1]-Y_range[0])*(Y.size()-1);
                if(index == (Y.size()-1)){
                    return index - 1;
                }
                else{
                    return index;
                } 
            }
            std::vector<Real> _a ;
            std::string X_name, Y_name , f_name , Filebase;
            double X_range[2] , Y_range[2];
    
    };

    /// This file only includes the functions requred for reading values from a set of coefficients. To set up the coefficients, 2Dinterpolator.cpp is required.

    template<class Real>
    void BicubicInterpolator<Real>::setup_from_coeff(std::string &filebase,std::string &x_name , std::string &y_name , std::string &F_name , int *grid_size){
        X_name = x_name;
        Y_name = y_name;
        f_name = F_name;
        read_XY(filebase, grid_size);

        X_range[0] = *std::min_element(X.begin(),X.end());
        X_range[1] = *std::max_element(X.begin(),X.end());
        Y_range[0] = *std::min_element(Y.begin(),Y.end());
        Y_range[1] = *std::max_element(Y.begin(),Y.end());
        read_coeff(filebase);
        Filebase = filebase;
    }

    template<class Real>
    void BicubicInterpolator<Real>::read_coeff(std::string &filebase){
        std::ifstream fin;
        std::cout << "Reading coefficients..." << std::endl;
        std::string coeffile = filebase+"/coeff/"+f_name+std::to_string(X.size())+"x"+std::to_string(Y.size())+".txt";
        
        fin.open(coeffile);
            if(!fin.is_open()){
                file_error(coeffile);
            }
            std::cout << coeffile << std::endl;
            Real line;
            while(fin >> line){
		        _a.push_back(line);
	        }
        fin.close();
        std::cout << "Coefficients read " << std::endl;
    }

    template<class Real>
    void BicubicInterpolator<Real>::read_XY(std::string &filebase, int *grid_size){
        
        std::ifstream fin;
        Real line;
        std::string gridfile = filebase+"/grid/"+X_name+std::to_string(grid_size[0])+".txt";
        fin.open(gridfile);
        if(!fin.is_open()){
              std::string error_message("\nCoefficient file not found or file error!\nFilename: \n ");
            error_message.append( gridfile ) ;
            throw std::invalid_argument(error_message );
        }
            while(fin >> line){
		        X.push_back(line);
	        }
        fin.close();
        gridfile = filebase+"/grid/"+Y_name+std::to_string(grid_size[1])+".txt"; 
        fin.open(gridfile);
        if(!fin.is_open()){
            std::string error_message("\nCoefficient file not found or file error!\nFilename: \n ");
            error_message.append( gridfile ) ;
            throw std::invalid_argument(error_message );
        }   
            while(fin >> line){
		        Y.push_back(line);
	        }
        fin.close();

    }
  
    template<class Real>
    int*
    BicubicInterpolator<Real>::return_index(Real x ,Real y)
    {
        // no extrapolation...
        if( x < *std::min_element(X.begin(),X.end())
        || x > *std::max_element(X.begin(),X.end())
        || y < *std::min_element(Y.begin(),Y.end())
        || y > *std::max_element(Y.begin(),Y.end())
        ){
            std::string error_message = "\nValue out of range of interpolator: ";
            error_message.append(Filebase);
            char out_of_range_values[150];
            snprintf(out_of_range_values, 150 , "\nx = %.3g, y = %.3g\nxmin = %.3g , xmax = %.3g\nymin = %.3g , ymax = %.3g" , x,y , X_range[0] ,X_range[1] , Y_range[0], Y_range[1]);
            error_message = error_message+out_of_range_values;
            throw std::invalid_argument(error_message);
            return 0;
        }

        static int indexes[2];
        indexes[0] = this->get_x_index_to_left_of(x);
        indexes[1] = this->get_y_index_below(y);
        
        if(indexes[0] < 0)
            indexes[0] = 0;
        if(indexes[1] < 0)
            indexes[1] = 0;
  

        return indexes;
    }

    template<class Real>
    Real
    BicubicInterpolator<Real>::return_value(Real x ,Real y , int i , int j)
    {
        Real *a_data = &_a.front();
        Eigen::Map<Matrix44>a(a_data + 16*(Y.size()-1)*i + 16*j, 4, 4) ; 
       
        Real xL = X[i+1] - X[i];
        Real yL = Y[j+1] - Y[j];

        // now, create the coordinate vectors (see Wikipedia https://en.wikipedia.org/wiki/Bicubic_interpolation)
        RowVector4 vx;
        vx[0] = 1;                                   // x^0
        vx[1] = (x - X[i])/xL;                    // x^1
        vx[2] = vx[1] * vx[1];                       // x^2
        vx[3] = vx[2] * vx[1];                       // x^3

        ColVector4 vy;
        vy[0] = 1;                                   // y^0
        vy[1] = (y - Y[j])/yL;                    // y^1
        vy[2] = vy[1] * vy[1];                       // y^2
        vy[3] = vy[2] * vy[1];                       // y^3


        // interpolation is just x*a*y

        return vx*a*vy;
    }

    template<class Real>
    Real*  BicubicInterpolator<Real>::return_grads(Real x ,Real y , int i , int j )
    {
        static Real grads[2];

        Real *a_data = &_a.front();
        Eigen::Map<Matrix44>a(a_data + 16*(Y.size()-1)*i + 16*j, 4, 4) ; 
      
        Real xL = X[i+1] - X[i];
        Real yL = Y[j+1] - Y[j];

        // now, create the coordinate vectors (see Wikipedia https://en.wikipedia.org/wiki/Bicubic_interpolation)
        
        //x gradient
        RowVector4 vx;
        ColVector4 vy;

        vx[0] = 0;                                  //0
        vx[1] = 1.0/xL;                                   // x^0
        vx[2] = 2*(x - X[i])/(xL*xL);                    // x^1
        vx[3] = 3.0/4*vx[2] * vx[2] *xL;                       // x^2
        
        vy[0] = 1;                                   // y^0
        vy[1] = (y - Y[j])/yL;                    // y^1
        vy[2] = vy[1] * vy[1];                       // y^2
        vy[3] = vy[2] * vy[1];                       // y^3

        grads[0] = vx*a*vy;

        //y gradient

        vx[0] = 1;                                   // x^0
        vx[1] = (x - X[i])/xL;                    // x^1
        vx[2] = vx[1] * vx[1];                       // x^2
        vx[3] = vx[2] * vx[1];                       // x^3

       
        vy[0] = 0;                                  // 0 
        vy[1] = 1.0/yL;                                   // y^0
        vy[2] = 2*(y - Y[j])/(yL*yL);                    // y^1
        vy[3] = 3.0*((y - Y[j])/yL)*((y - Y[j])/yL)/yL;                       // y^2
        

        grads[1] = vx*a*vy;

        return grads;
        }

}

#endif