#ifndef _Core_Latent_H
#define _Core_Latent_H

#include <Eigen/Dense>

/// @brief Class for calculating the latent heat due to solidification of the core
class Core_Latent{
    public:
        Eigen::ArrayXd Heat_old(double *solid_mass , Eigen::ArrayXd core_is_liquid , double &dt){
            
            Eigen::ArrayXd Latent_heat = Eigen::ArrayXd::Zero(core_is_liquid.size());
            
            if(abs(solid_mass[1]-solid_mass[0]) > 0.0){ 
                Latent_heat = DeltaL*(1.0-core_is_liquid)*(solid_mass[1]-solid_mass[0])/solid_mass[1]/dt;
                std::cout << "calcing latent " << solid_mass[0] << " solids " << solid_mass[1] << std::endl;
            }               
            return Latent_heat;
        }

        Eigen::ArrayXd Heat(double *solid_mass , const double &core_mass , int &Jcore , double &dt){
            Eigen::ArrayXd Latent_heat = Eigen::ArrayXd::Zero(Jcore);
            if(abs(solid_mass[1]-solid_mass[0]) > 0.0){ 
                Latent_heat = Eigen::ArrayXd::Constant(Jcore , DeltaL*(solid_mass[1]-solid_mass[0])/core_mass/dt);
                std::cout << "calcing latent " << solid_mass[0] << " solids " << solid_mass[1] << std::endl;
            }               
            return Latent_heat;
        }
        
    private:
        double DeltaL = 247e3;
        
};

#endif