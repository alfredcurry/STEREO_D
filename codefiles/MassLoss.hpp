#ifndef _MASS_LOSS_H
#define _MASS_LOSS_H

#include "Vector_reader.h"
#include "2Dinterpolator.hpp"
#include "Kang_analytical.hpp"

template <typename Real>
Real Fortney_R(Real const &logM , Real const &RMF){
    //Note CMF of 0.9 gives rho_fac ~0.6 0.8 ~0.65 0.7 ~ 0.7 0.99 ~ 0.58 (cold)
    return (0.0592*RMF + 0.0975)*pow(logM,2) + (0.2337*RMF + 0.4938)*logM + (0.3102*RMF + 0.7932);
}


template <typename Real>
class MassLoss
{
    public:
        static Real Constant(){
            return 1e7;
        }
        static Real Isothermal(Real const &M){
            // Approx to isothermal solution in Perez-Becker and Chiang (Fig.2). Takes input in M_Earth, gives result in M_\Earth / Gyr
            double a , M_crit = 3e-3 , M_dot_peak = 20;
            Real A, B;
         
            a = log(M_dot_peak/1.0e-1)/7e-2;
            A = M_dot_peak*exp(a*M_crit);
            B = M_dot_peak/pow(M_crit,2.0/3);

            if( M > M_crit ){
                return A*exp(-a*M);
            }else{
                return B*pow(M,2.0/3);
            }
        }
        static Real Perez_PowerLaw(Real const &M){
            // Approx to isothermal solution in Perez-Becker and Chiang (Fig.2). Takes input in M_Earth, gives result in M_\Earth / Gyr
            double M_crit = 3e-3 , M_dot_peak, a = -3.5;
            Real A, B;

            A = 1.8e-4;         
            M_dot_peak = A*pow(M_crit,a);
            B = M_dot_peak/pow(M_crit,2.0/3);

            if( M > M_crit ){
                return A*pow(M,a);
            }else{
                return B*pow(M,2.0/3);
            }
        }
        Real PBC_tabulated(Real const &M){
            // Tabulated solution form Perez-Becker and Chiang (Fig.8). Takes input in M_Earth, gives result in M_\Earth / Gyr
            return pow(10 , Mass_Loss_spline(log10(M)));
        }
        boost::math::interpolators::cardinal_cubic_hermite<std::vector<Real>> Mass_Loss_spline;
        
        Real Bicubic_tabulated(Real const &M , Real const &R , Real &rho_fac){
            // Tabulated solution from Booth (2022). Takes input in M_Earth and R_Earth gives result in kg/s
            int *indexes;
            
            indexes = Mass_Loss_Bicubic.return_index(log10(M) , rho_fac);  

            return pow(10.0 , Mass_Loss_Bicubic.return_value(log10(M) , rho_fac ,indexes[0],indexes[1]))/1e3;
        }
        
        Real Kang_anal(Real const &g , Real const &R , Real const &period ){
            return K_mdot.day_side(g,R,period) + K_mdot.night_side(g,R);
        }
        Real Kang_terminator(Real const &g , Real const &R ){
            return K_mdot.terminator(g,R);
        }
        void read_table(char &type , double &Tss){
            std::string file , filebase = "MassLossTables/";
            
            std::string comp = "SiO";
            K_mdot.set_composition(comp);
            K_mdot.set_Tpool(Tss);

            if(type == 'p'){
                std::cout << "PBC mass loss" << std::endl;

                int ROWS = 64;
                double M0 = -2.2, dM = 0.02;
                std::vector<Real> MassLoss_data(ROWS) , MassLoss_dM_data(ROWS) ;
                char filestring[100];
                snprintf(filestring , 100 , "Perez-Becker_Chiang/M_dot%.d_%.2g_%.1g.txt" , ROWS , M0 , dM);
                std::cout << "Reading Mass Loss..." << std::endl;
                file = filebase+filestring;
                readVector(MassLoss_data , ROWS , file );
                
                snprintf(filestring , 100 , "Perez-Becker_Chiang/M_dotdM%.d_%.2g_%.1g.txt" , ROWS , M0 , dM);
                std::cout << "Reading Mass Loss derivative..." << std::endl;
                file = filebase+filestring;
                readVector(MassLoss_dM_data , ROWS , file);

                Mass_Loss_spline = boost::math::interpolators::cardinal_cubic_hermite<std::vector<Real>>(std::move(MassLoss_data), std::move(MassLoss_dM_data), M0, dM);

            }else if( type == 'B'){
                std::cout << "Booth mass loss" << std::endl;
                filebase += "Booth/";
                if(Tss == 1460){
                    filebase += "a0.0337Tss1463.8";
                }else if(Tss == 1640){
                    filebase += "a0.0267Tss1642.4";
                }else if(Tss == 1840){
                    filebase += "a0.0212Tss1842.8";
                }else if(Tss == 1950){
                    filebase += "a0.0189Tss1952.0";
                }else if(Tss == 2070){
                    filebase += "a0.0169Tss2067.7";
                }else if(Tss == 2190){
                    filebase += "a0.0150Tss2190.2";
                }else if(Tss == 2320){
                    filebase += "a0.0134Tss2320.0";
                }else if(Tss == 2460){
                    filebase += "a0.0119Tss2457.5";
                }else if(Tss == 2600){
                    filebase += "a0.0106Tss2603.1";
                }else{
                    throw std::invalid_argument("No Mass Loss for this Tss");
                }
                std::string x_name , y_name ;
                file = "logmdot";
                int grid_size[2];
                grid_size[0] = 20;
                grid_size[1] = 7;
                x_name = "log_m";
                y_name = "rho_fac";
                Mass_Loss_Bicubic.setup_from_coeff(filebase, x_name, y_name, file, grid_size);
                
            }
        }
    private:
        KangMdot<Real> K_mdot;
        Interpolators::BicubicInterpolator<double> Mass_Loss_Bicubic;

};

#endif