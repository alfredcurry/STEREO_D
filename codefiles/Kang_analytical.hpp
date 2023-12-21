#ifndef _KANG_H
#define _KANG_H

#include <cmath>
#include <stdexcept>
#include <boost/math/tools/roots.hpp>
#include "PhysicalConsts.h"

template <typename Real> /// Class for calculting mass loss using the analytical model from Kang (2021)
class KangMdot
{
    public:
        Real day_side( Real const &g , Real const &Rp , Real const &period ){
            Real T_sonic = T_sonic_func(g , Rp , period);
            return (M_PI * A_sat * g*g * pow(Rp,4) * (1-cos_theta_pool) ) /(2 * pow(R_gas_spec * T_sonic , 2.5)) * exp(-B_sat/T_sonic);
        }
        Real terminator( Real const &g , Real const &Rp ){
            return 2*M_PI * P_chem(T_pool) * sqrt(R_gas_spec*T_pool) * sin_theta_pool * Rp/g ;
        }
        Real antistellar_terminator( Real const &g , Real const &Rp){
            return exp(-M_PI*g*Rp / (4*sqrt(Cp*T_pool * M_PI * R_gas_spec*T_sat)));
        }
        Real night_side_antistellar( Real const &g, Real const &Rp){
            return exp(-g*Rp/(latent_heat/mol_mass));
        }
        Real night_side( Real const &g , Real const &Rp ){
            return terminator(g , Rp )*antistellar_terminator(g , Rp)*night_side_antistellar(g, Rp);
        }
        void set_Tpool(Real &T_ss , std::string const &surface_profile="cos025"){
            Real T_c = 1673;
            calculate_T_pool(T_ss , T_c );
            T_sat = calculate_T_sat();
            std::cout << T_sat << " cos " << cos_theta_pool << " sin " << sin_theta_pool << std::endl;
        }
        void set_composition(std::string &composition){
            if (composition == "Na"){
                mol_mass = 0.023 ; //JK-1kg-1
                latent_heat = 96.96e3 ; //Jmol-1
                Cp = 903.3 ; //JK-1kg-1
                A_chem = pow(10,9.6) ; //Pa
                B_chem = 38000 ; //K
                A_sat = pow(10,9.54) ; //Pa
                B_sat = 12070.4 ; //K
            }else
            if (composition == "SiO"){
                mol_mass = 0.044; //kgmol-1
                latent_heat = 411.5e3 ; //Jmol-1
                Cp = 661 ; //JK-1kg-1
                A_chem = pow(10,14.086) ; //Pa
                B_chem = 70300 ; //K
                A_sat = pow(10,13.1) ; //Pa
                B_sat = 49520 ; //K
            }else
            if ( composition == "Al2O3"){
                mol_mass = 0.102 ; //kgmol-1
                latent_heat = 643e3 ; //Jmol-1
                Cp = 755 ; //JK-1kg-1
                A_chem = pow(10,12.5) ; //Pa
                B_chem = 71380 ; //K
                A_sat = pow(10,16.1) ; //Pa
                B_sat = 77365 ; //K
            }else{
                throw std::invalid_argument("Invalid flow composition.");
            }
            R_gas_spec = R_gas / mol_mass ; //JK-1kg-1
            adiabatic_kappa = R_gas_spec/Cp;
        }
    private:
        Real R_gas_spec , mol_mass , latent_heat , Cp , A_chem , B_chem , A_sat , B_sat , adiabatic_kappa , T_sat , cos_theta_pool , sin_theta_pool , T_pool ;


        Real T_sonic_func( Real const &g , Real const &Rp , Real const &period ){
            return T_sat * exp( (grav_potential(Rp,g,Rp,period) - grav_potential(roche_R(Rp,g,period),g,Rp,period) + Cp * (T_pool - T_sat))/(R_gas_spec*B_sat) );
        }
        Real grav_potential(Real const &r , Real const &g , Real const &Rp , Real const &period){
            return -g*Rp*Rp/r - 1.5 * pow(2*M_PI/period,2) * r*r;
        }
        Real roche_R(Real const &Rp , Real const &g , Real const &period){
            return pow(g*Rp*Rp / (3*pow(2*M_PI/period,2)),1.0/3);
        }
        Real P_chem( Real &T ){
            return P_chemical_general(A_chem,B_chem)(T);
        }
        Real P_sat( Real &T ){
            return P_chemical_general(A_sat,B_sat)(T);
        }
        struct P_chemical_general
        {
            P_chemical_general(Real const &_A , Real const &_B):
                A(_A) , B(_B)
            {}
            Real operator()(Real const &T){
                return A * exp(-B/T);
            }
            private:
                Real A , B;
        };
        struct T_sat_functor
        {
            T_sat_functor(Real const &_adiabatic_kappa , Real const &_T_pool , Real const &_A_sat ,Real const &_B_sat , Real const &_A_chem ,Real const &_B_chem):
                adiabatic_kappa(_adiabatic_kappa) , T_pool(_T_pool) , A_sat(_A_sat) , B_sat(_B_sat) , A_chem(_A_chem) , B_chem(_B_chem)
            {}
            Real operator()(Real const &T_sat){    
                return T_sat * pow(P_chemical_general(A_sat,B_sat)(T_sat),-adiabatic_kappa) - T_pool * pow(P_chemical_general(A_chem,B_chem)(T_pool),-adiabatic_kappa);
            }
            private:
                Real adiabatic_kappa , T_pool , A_sat , B_sat , A_chem , B_chem;
        };   
        Real calculate_T_sat(){
             // Solve using bracket_and_solve (no derivatives).
            using namespace boost::math::tools;           // For bracket_and_solve_root.

            const boost::uintmax_t maxit = 50;            // Limit to maximum iterations.
            boost::uintmax_t it = maxit;                  // Initially our chosen max iterations, but updated with actual.
            int digits = std::numeric_limits<Real>::digits;  // Maximum possible binary digits accuracy for type T.
            // Some fraction of digits is used to control how accurate to try to make the result.
            int get_digits = digits - 0;      
                         
            eps_tolerance<Real> tol(get_digits);             // Set the tolerance.

            Real lower_bound = T_pool/10.0;
            std::cout << T_pool  << " T pool " << adiabatic_kappa << std::endl;
            std::pair<Real, Real> r = toms748_solve(T_sat_functor(adiabatic_kappa , T_pool, A_sat , B_sat , A_chem , B_chem), lower_bound, T_pool , tol, it);
  
            if (it >= maxit)
            {
                std::cout << "Unable to locate solution with Brent's method in " << maxit << " iterations:"
                " Current best guess is between " << r.first << " and " << r.second << std::endl;
            }

            return r.first + (r.second - r.first)/2;      // Midway between brackets is our result, if necessary we could
                                                // return the result as an interval here.                                               
        }
        void calculate_T_pool( Real &T_ss , Real &T_c , std::string const &surface_profile="cos025"){
            if( surface_profile == "cos025"){
                cos_theta_pool = pow(T_c/T_ss,4);
                if( cos_theta_pool > 1.0){
                    throw std::invalid_argument("Substellar temperature less than melting point.");
                }
                sin_theta_pool = sqrt(1.0-cos_theta_pool*cos_theta_pool);

                T_pool = 4.0/5 * T_ss * (1 - pow(cos_theta_pool,5.0/4)) / (1 - cos_theta_pool);
            }else{
                throw std::invalid_argument("Unknown surface temperature profile.");
            }
        }
};


#endif