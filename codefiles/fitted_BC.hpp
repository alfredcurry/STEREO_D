#include <iostream>
#include "PhysicalConsts.h"
#include <vector>

/// @brief Class containing the functional form of the fitted function for the asymmetric boundary condition used for highly heated, tidally locked planets
class Fitted_Layer{

    public:
        double operator()( const double &F , const double &g){/// The function
            return pow(pow( pow( pow(power_law(F , g , A1 ,a1 , b1) + Tend , -beta) + pow( log_law(F, g ,A2 , a2 , b2), -beta ) , gam/beta) + pow(Tcr , -gam) ,-alpha/gam) + pow(power_law(F , g , A4 , a4 , b4) , alpha ) , 1.0/alpha);
        }
        double g_derivative(const double &F , const double &g){/// derivative w.r.t. g (the surface gravity)
            return pow(pow( pow( pow(power_law(F , g , A1 ,a1 , b1) + Tend , -beta) + pow( log_law(F, g ,A2 , a2 , b2), -beta ) , gam/beta) + pow(Tcr , -gam) ,-alpha/gam) + pow(power_law(F , g , A4 , a4 , b4) , alpha ) , 1.0/alpha - 1) * (
                pow(power_law(F , g , A4 , a4 , b4), alpha -1)*power_lawdg(F , g , A4 , a4 , b4) 
                    + pow( pow( pow(power_law(F , g , A1 ,a1 , b1) + Tend , -beta) + pow( log_law(F, g ,A2 , a2 , b2), -beta ) , gam/beta) + pow(Tcr , -gam) ,-alpha/gam - 1) 
                        * pow( pow(power_law(F , g , A1 ,a1 , b1) + Tend , -beta) + pow( log_law(F, g ,A2 , a2 , b2), -beta ) , gam/beta -1) 
                            * (pow(power_law(F , g , A1 ,a1 , b1) + Tend , -beta-1)*power_lawdg(F , g , A1 ,a1 , b1) + pow( log_law(F, g ,A2 , a2 , b2), -beta -1 ) *log_lawdg(g,b2) )
            );
        }
        double F_derivative(const double &F , const double &g){/// derivative w.r.t. F (the mean flux through P0 / to the surface)
            return pow(pow( pow( pow(power_law(F , g , A1 ,a1 , b1) + Tend , -beta) + pow( log_law(F, g ,A2 , a2 , b2), -beta ) , gam/beta) + pow(Tcr , -gam) ,-alpha/gam) + pow(power_law(F , g , A4 , a4 , b4) , alpha ) , 1.0/alpha - 1) * (
                pow(power_law(F , g , A4 , a4 , b4), alpha -1)*power_lawdF(F , g , A4 , a4 , b4) 
                    + pow( pow( pow(power_law(F , g , A1 ,a1 , b1)+ Tend , -beta) + pow( log_law(F, g ,A2 , a2 , b2), -beta ) , gam/beta) + pow(Tcr , -gam) ,-alpha/gam - 1) 
                        * pow( pow(power_law(F , g , A1 ,a1 , b1)+ Tend , -beta) + pow( log_law(F, g ,A2 , a2 , b2), -beta ) , gam/beta -1) 
                            * (pow(power_law(F , g , A1 ,a1 , b1)+ Tend , -beta-1)*power_lawdF(F , g , A1 ,a1 , b1) + pow( log_law(F, g ,A2 , a2 , b2), -beta -1 ) *log_lawdF(F,a2) )
            );
        }
        void Assign(std::vector<double> &_variables){
            A1 = _variables[0];
            a1 = _variables[1];
            b1 = _variables[2]; 
            A2 = _variables[3]; 
            a2 = _variables[4]; 
            b2 = _variables[5]; 
            A4 = _variables[6]; 
            a4 = _variables[7]; 
            b4 = _variables[8]; 
            beta = _variables[9]; 
            alpha = _variables[10]; 
            gam = _variables[11]; 
            Tend = _variables[12];
            Tcr = _variables[13];
        }
    private:
        double A1 , a1 , b1 , A2 , a2 , b2 , A4 , a4 , b4 , beta , alpha , gam , Tcr , Tend ;
        double Tnight(const double &F){
            return pow( F/Sig_Stefan , 0.25);
        }
        double TnightdF(const double &F){
            return 0.25*pow( F , -0.75)/pow(Sig_Stefan , 0.25);
        }
        double power_law(const double &F , const double &g , double &A , double &a , double &b){
            return Tnight(F) + A*pow(F,a) * pow(g,b);
        }
        double power_lawdF( const double &F , const double &g , double &A , double &a , double &b){
            return TnightdF(F) + a*A*pow(F,a-1) * pow(g,b);
        }
        double power_lawdg( const double &F , const double &g , double &A , double &a , double &b){
            return  b*A*pow(F,a) * pow(g,b-1);
        }
        double log_law(const double &F , const double &g , double &A , double &a , double &b){
            return A + log(F)*a + log(g)*b ;
        }
        double log_lawdF(const double &F , const double &a ){
            return a/F;
        }
        double log_lawdg(const double &g , const double &b){
            return b/g ;
        }
};  