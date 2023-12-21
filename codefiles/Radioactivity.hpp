#ifndef _Radioactivity_H
#define _Radioactivity_H

#include <Eigen/Dense>

/// @brief Class for calculating the power per unit mass from radioactivity
class Radioactivity{
    public:
        Eigen::ArrayXd Constant(int J){
            return Eigen::ArrayXd::Constant(J,H);
        }
        Eigen::ArrayXd Earth(int J , double t){
            H = H_Th*exp(-(t-t0)/t_Th) + H_U238*exp(-(t-t0)/t_U238) + H_K*exp(-(t-t0)/t_K) + H_U235*exp(-(t-t0)/t_U235);
            
            return Eigen::ArrayXd::Constant(J,H);
        }
        Eigen::ArrayXd Iron60_Earth(int J , double t){
            H = SLR(t,t_Fe,H_Fe) ;            
            return Eigen::ArrayXd::Constant(J,H);
        }
        Eigen::ArrayXd Al26_Earth(int J , double t){
            H = SLR(t,t_Al,H_Al) ;
            return Eigen::ArrayXd::Constant(J,H);
        }
        double SLR(double t , double t_SLR , double H_SLR){
            return H_SLR*exp(-t/t_SLR) ;
        }
        void set_Constant(double &_H){
            H = _H;
        }
    private:
        double H;
        double H_Th = 3.27e-12 , H_U238 = 2.91e-12 , H_K = 1.08e-12 , H_U235 = 0.125e-12; //Wkg-1
        double t_Th = 4.47e9/log(2) , t_U238 = 4.47e9/log(2) , t_K = 1.25e9/log(2) , t_U235 = 0.704e9/log(2), t0 = 4.5e9; //yrs
        double H_Fe = 6.28e-11 , t_Fe = 2.62e6/log(2); //Lugaro (2018)
        double H_Al = 1.55e-7 , t_Al = 1.035e6;//Lugaro (2018)
        
};

#endif