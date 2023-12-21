#include <iostream>
#include <Eigen/Dense>
#include "converge.hpp"
#include "StructStruct.hpp"
#include "PhysicalConsts.h"
#include "EquationOfState.hpp"
#include "MeltFraction.hpp"
#include "Radioactivity.hpp"
#include <boost/math/interpolators/pchip.hpp>
#include <vector> 
#include "MassLoss.hpp"
#include "MassTerm.hpp"
#include "Edge_layer.hpp"
#include "iron_core.hpp"
#include "Conductivity.hpp"
#include "Print_and_save.hpp"
#include "Interpolator.hpp"
#include "Core_Latent.hpp"

using namespace Eigen;
class stepper{
    public:
        Saver print;
        void k_alter(struct P_P *planet , EoS *eos , double k_aim, double acc , char *num  , std::string &_folder){
            /** Function for altering the conductivity, k, in small increments */
            con.nmax = 50;
            print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit);
            static_struct(planet , eos , acc);
            print.save_values(*num , planet , eos  , time , f_time , con);
            double k_step;
            con.nmax = 100;
            while (abs(planet->conds[0] - k_aim) != 0 ){
                k_step = std::min(abs(k_aim* 0.05), abs(k_aim-planet->conds[0]));
                if( k_aim < planet->conds[0]){k_step = -k_step;} 
                planet->conds[0] += k_step;
                std::cout << k_aim << " / " << planet->conds[0] << std::endl;
                static_struct(planet, eos , acc);
            }
            print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit);
            print.save_values(*num , planet , eos  , time , f_time , con);

        }
	    void n_alter(struct P_P *planet , EoS *eos , double n_aim, double acc , char *num  , std::string &_folder){
            /** Function for altering the adiabatic index, n,  in small increments */
            con.nmax = 50;
            print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit);
            static_struct(planet , eos , acc);
            print.save_values(*num , planet , eos  , time , f_time , con);
            double n_step;
            con.nmax = 100;
            while (abs(planet->n - n_aim) != 0 ){
                n_step = std::min(abs(n_aim* 0.005), abs(n_aim-planet->n));
                if( n_aim < planet->n){n_step = -n_step;} 
                planet->n += n_step;
                std::cout << n_aim << " / " << planet->n << std::endl;
                static_struct(planet , eos , acc);
            }
            print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit);
		print.save_values(*num , planet , eos  , time , f_time , con);

	    }
        void tau_alter(struct P_P *planet , EoS *eos , double tau_aim, double acc , char *num  , std::string &_folder){
            /** Function for altering the atmsopsphere optical depth, tau, in small increments */
            con.nmax = 200;
            print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit);
            con.converge(planet , eos , acc);
            double tau_step;
            con.nmax = 100;
            while (abs(planet->tau- tau_aim) != 0 ){
                tau_step = std::min(abs(tau_aim * 0.00001), abs(tau_aim-planet->tau));
                if( tau_aim < planet->tau){tau_step = -tau_step;} 
                planet->tau += tau_step;
                std::cout << "tau: " << tau_aim << " / " << planet->tau << std::endl;
                con.converge(planet , eos , acc);
                print.print_accuracies(planet , t , time , con );

            if(con.maxinacc > acc){
                print.Not_Converged(planet , eos , con);            }
            }
            print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit);
            print.save_values(*num , planet , eos  , time , f_time , con);

        }
        void L_alter(struct P_P *planet , EoS *eos , double L_aim, double acc , char *num  , std::string &_folder){
            /** Function for altering the planet luminosity, L, in small increments */

            con.nmax = 1000;
            print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit);
            std::cout << "L "  << planet->scale(3) << " " << planet->L_core << std::endl;    
            char nabtype = 'z';
            if (planet->nabtype == 'c' || planet->nabtype == 'g' || planet->nabtype == 'B'){
                nabtype = planet->nabtype;
                planet->nabtype = 'a';
            }
	        static_struct(planet , eos , acc);
            print.save_values(*num , planet , eos  , time , f_time , con);
            double L_step;

            while (abs(planet->scale(3) - L_aim) != 0 ){
                L_step = std::min(abs(planet->scale(3) * 0.001), abs(L_aim-planet->scale(3)));
              
		    if( L_aim < planet->scale(3)){L_step = -L_step;} 
                planet->scale(3) += L_step;
                planet->L_core = planet->scale(3);
                static_struct(planet , eos , acc);
            }
            print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit);
            print.save_values(*num , planet , eos  , time , f_time , con);
            planet->nabtype = nabtype ;

        }
        /**
         * I can't remember the difference for this one compared to above
        void L_core_cheat(struct P_P *planet , EoS *eos , double L_aim, double acc , char *num  , std::string &_folder){
            con.nmax = 1000;
            print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit);
            std::cout << "L "  << planet->scale(3) << " " << planet->L_core << std::endl;    
	        static_struct(planet , eos , acc);
            print.save_values(*num , planet , eos  , time , f_time , con);
            double L_step;
       
            while (abs(planet->scale(3) - L_aim) != 0 ){
                L_step = std::min(abs(planet->scale(3) * 0.01), abs(L_aim-planet->scale(3)));
              
		    if( L_aim < planet->scale(3)){L_step = -L_step;} 
                planet->scale(3) += L_step;
                planet->L_core = planet->scale(3);
                static_struct(planet , eos , acc);
            }
            print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit);
            print.save_values(*num , planet , eos  , time , f_time , con);

        }*/
        void temp_alter(struct P_P *planet, EoS *eos , double T_aim, double acc , char *num  , std::string &_folder){
            /** Function for altering the planet outer temperature, Tsurf, in small increments */

            f_time = 0;
            con.nmax = 150;
            print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit);
            static_struct(planet , eos, acc);
            print.save_values(*num , planet , eos  , time , f_time , con);
            double T_step;
            double Tsurfnext , print_T;
            if(planet->Tsurf[0] > T_aim){
                print_T = -50;
            }else{
                print_T = 50;
            }
            Tsurfnext = planet->Tsurf[0] + print_T ;
            while (abs(planet->Tsurf[0] - T_aim) != 0 ){
                T_step = std::min(abs(planet->Tsurf[0] * 1e-3), abs(T_aim-planet->Tsurf[0]));
                if( T_aim < planet->Tsurf[0]){T_step = -T_step;} 
                planet->Tsurf[0] += T_step;

                if( (planet->Tsurf[0] > Tsurfnext && print_T > 0) || (planet->Tsurf[0] < Tsurfnext && print_T < 0 )){
                    planet->Tsurf[0] = Tsurfnext;
                }         
                print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit);

                static_struct(planet , eos , acc);

                ArrayXXd problems(J,2);
                problems.col(0) = eos->melt.NormT;
                problems.col(1) = planet->phi;
                //std::cout << problems << std::endl;
            
                if(planet->Tsurf[0] == Tsurfnext){
                    print.save_values(*num , planet , eos  , time , f_time , con);
                    Tsurfnext = Tsurfnext+print_T ; 
                    if((Tsurfnext > T_aim && print_T > 0) || (Tsurfnext < T_aim && print_T < 0 )){
                        Tsurfnext = T_aim ;
                    }
                }
                
            }

        }
        void Teq_alter(struct P_P *planet, EoS *eos , double T_aim, double acc , char *num  , std::string &_folder){
            /** Function for altering the planet equilibrium temperature, Teq, in small increments */

            f_time = 0;
            con.nmax = 2000;
            print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit);
            static_struct(planet , eos , acc);
            print.save_values(*num , planet , eos  , time , f_time , con);
            double T_step;
            
            while (abs(planet->Teq - T_aim) != 0 ){
                T_step = std::min(abs(1.0), abs(T_aim-planet->Teq));
                
                if( T_aim < planet->Teq){T_step = -T_step;} 
                planet->Teq += T_step;
                static_struct(planet , eos , acc);
            }
            print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit);
            print.save_values(*num , planet , eos  , time , f_time , con);

        }
        void mu_m_alter(struct P_P *planet, EoS *eos , double mu_aim, double acc , char *num  , std::string &_folder){
            /** Function for altering the mean molecuar weight of ideal gas planet, mu, in small increments */

            f_time = 0;
            con.nmax = 500;
            print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit);
            static_struct(planet , eos , acc);
            print.save_values(*num , planet , eos  , time , f_time , con);
            
            double mu_step;
   
            while (abs(planet->mu_m  - mu_aim) != 0 ){
                mu_step = std::min(abs(planet->mu_m* 0.1), abs(mu_aim-planet->mu_m));
                if( mu_aim < planet->mu_m){mu_step = -mu_step;} 
                planet->mu_m += mu_step;
                static_struct(planet , eos , acc);
            }
            print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit);
            print.save_values(*num , planet , eos  , time , f_time , con);

        }
        void pressure_alter(struct P_P *planet, EoS *eos , double P_aim, double acc , char *num  , std::string &_folder){
            /** Function for altering the planet outer pressure, Psurf, in small increments */

            f_time = 0;
            con.nmax = 500;
            print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit);
            static_struct(planet , eos , acc);
            print.save_values(*num , planet , eos  , time , f_time , con);
            planet->Pbound = 'f';
            planet->Tbound = 'f';
            double P_step;
   
            while (abs(planet->Psurf[0] - P_aim) != 0 ){
                P_step = std::min(abs(planet->Psurf[0] * 0.01), abs(P_aim-planet->Psurf[0]));
                if( P_aim < planet->Psurf[0]){P_step = -P_step;} 
                planet->Psurf[0] += P_step;
                static_struct(planet , eos , acc);
            }
            print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit);
            print.save_values(*num , planet , eos  , time , f_time , con);

        }
        void simp_M_alter(struct P_P *planet, EoS *eos , double M_aim, double M_interval , double acc , char *num  , std::string &_folder){
            /** Function for altering the total planet mass in small increments (not time evolving)*/

            print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit);
            con.nmax = 500;
            static_struct(planet , eos , acc);
            print.save_values(*num , planet , eos  , time , f_time , con);
            std::cout << " M aim " << M_aim/M_Earth << std::endl;
            double M_step , M_next;
            bool save;
            if(M_aim < planet->M){
                M_next = planet->M - M_interval*M_Earth;
            }else{
                M_next = planet->M + M_interval*M_Earth;

            }
            while(abs(planet->M - M_aim)/M_Earth > 1e-10){
                M_step= planet->M*0.1 ;
                M_step = std::min(abs(M_step) , abs(M_aim-planet->M));
                if( M_aim < planet->M){M_step = -M_step;}
                planet->M +=M_step;
                if(M_step < 0){
                    if(planet->M < M_next){
                        planet->M = M_next;
                        if (planet->M < M_aim){
                            planet->M = M_aim;
                        }else{
                            M_next = M_next - M_interval*M_Earth;
                        }               
                        save = true;
                    }
                }else{
                    if(planet->M > M_next){
                        planet->M = M_next;
                        if (planet->M > M_aim){
                            planet->M = M_aim;
                        }else{
                            M_next = M_next + M_interval*M_Earth;
                        }                        
                        save = true;
                    }
                }
                if(planet->M ==  M_aim){
                    save = true;
                }
                planet->M_total = planet->M;
                print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit);
                static_struct(planet, eos , acc);
                if(save == true){                    
                    print.save_values(*num , planet , eos  , time , f_time , con);
                }
                save = false;
            }
            
        }

        void core_incrementer(struct P_P *planet, EoS *eos , double M_core_aim, double X_env, double acc , char *num  , std::string &_folder){
            /**Function for introducing a core that is NOT resolved by the code. Core defined by hard-coded M-R relation, and mass fraction of resolved part*/


            print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit);
            con.nmax = 500;

            static_struct(planet , eos , acc);
            
            double M_step = planet->M*0.01 , env_step = 0.01*planet->M;
    
            while (planet->M_core < M_core_aim){
                M_step = std::min(M_step , M_core_aim-planet->M_core);                
		        env_step = std::min(env_step , M_core_aim/(1-X_env) - planet->M);

		        planet->M_core += M_step;
                planet->M = planet->M_core/(1-X_env);//+env_step;

                planet->R_core = R_Earth * pow(planet->M_core/M_Earth , 0.25);
                
                planet->m = JacobMaker::power_law_m(J , (planet->M - planet->M_core)/planet->M , 2.3 );
           
                JacobMaker::prepm(planet);
                static_struct(planet, eos , acc);
                
            }
            print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit);
            print.save_values(*num , planet , eos  , time , f_time , con);
        }

        void core_incrementer2(struct P_P *planet, EoS *eos , double R_core_aim , double core_frac, double acc , char *num  , std::string &_folder){
            /** Function for introducing a core that is NOT resolved by the code. Core defined by core-radius and core-mass-fraction */

            print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit);
            planet->g_surf = G_Newt*planet->M/(planet->scale(1)*planet->scale(1));
            con.nmax = 500;
            double M_core_aim = core_frac*planet->M;
            char nabtype = 'z';
            if (planet->nabtype == 'c' || planet->nabtype == 'g' || planet->nabtype == 'B'){
                nabtype = planet->nabtype;
                planet->nabtype = 'a';
            }
            static_struct(planet , eos , acc);
            print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit);
            print.save_values(*num , planet , eos  , time , f_time , con);
            double M_step ;
    
            while (planet->M_core < M_core_aim){
                if(planet->M_core > 0){
                    M_step = planet->M_core*0.05;
                }else{
                    M_step = planet->M*0.02;
                }
                M_step = std::min(M_step , M_core_aim-planet->M_core);                

		        planet->M_core += M_step;

                planet->R_core = R_core_aim*planet->M_core/M_core_aim;
                planet->L_core = planet->H(0)*planet->M_core;
                planet->m = JacobMaker::power_law_m(J , (planet->M - planet->M_core)/planet->M , 1.5 );
                planet->g_surf = G_Newt*planet->M/(planet->scale(1)*planet->scale(1));

                JacobMaker::prepm(planet);
                static_struct(planet, eos , acc);
                
            }
            print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit);
            print.save_values(*num , planet , eos  , time , f_time , con);
            planet->nabtype = nabtype;
        }

        void introduce_core(struct P_P *planet, EoS *eos , double core_frac, double acc , char *num  , std::string &_folder){
            /** Function for introducing a core that IS resolved by the code */

            print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit);
            planet->g_surf = G_Newt*planet->M/(planet->scale(1)*planet->scale(1));
            con.nmax = 500;
            
            char nabtype = planet->nabtype;

            if (planet->nabtype == 'c' || planet->nabtype == 'g' || planet->nabtype == 'B'){

                planet->nabtype = 'a';
            }
            std::cout << " T bound now " << planet->Tbound << std::endl;
            static_struct(planet , eos , acc);
            print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit);
            print.save_values(*num , planet , eos  , time , f_time , con);
            
            int j=0;
            std::cout << "Intrducing core" << std::endl;
            while(planet->m(j) <= core_frac && j < J-1){
                planet->Jcore = j;
                j++;
            }

            planet->solid_core_mass[0] = 0.0;
            planet->solid_core_mass[1] = 0.0;
            planet->solid_core_mass[2] = 0.0;

            planet->M_core = planet->M*planet->m(planet->Jcore);


            ArrayXXd y_temp;
            ArrayXd m_temp , mm_temp;
            int J = planet->rho.size() , Jin = J/6;
            if( 1 == 1){
                std::cout << "regridding " << std::endl;
                y_temp = planet->y;
                m_temp = planet->m ;
                mm_temp = planet->mm;
                planet->m.tail(J-planet->Jcore) = JacobMaker::power_law_m(J-planet->Jcore,1.0-planet->m(planet->Jcore),1.0);
                j=0;
                //std::cout << planet->m.transpose() << std::endl;
              
              
                double gamma = 1.0 , dm0 = 5e-5 , dm = dm0;
                while (planet->m(planet->Jcore+j) + dm < planet->m(planet->Jcore+j+1))
                {   
                    j = j+1;
                    planet->m(planet->Jcore+j) = planet->m(planet->Jcore+j-1) + dm;
                    dm = dm0*pow(double(j),gamma);
                }
                std::cout << "got to " << j << std::endl;
                                
                JacobMaker::prepm(planet);

                planet->y.col(0) = interpolate(planet->mm , mm_temp , y_temp.col(0) ); // pressure, using interpolation
                planet->y.block(0,1,J-1,1) = interpolate(planet->m.tail(J-1) , m_temp , y_temp.col(1) ); //radius , using interpolation
                planet->y.col(2) = interpolate(planet->mm , mm_temp , y_temp.col(2) ); // pressure, using interpolation
                planet->y.block(0,3,J-1,1) = interpolate(planet->m.tail(J-1) , m_temp , y_temp.col(3) ); //radius , using interpolation
                
            ArrayXXd display(J,7);
            
            display.col(0) = ArrayXd::LinSpaced(J,0,J-1);
            display.col(2).tail(J-1) = planet->dm;
            display.col(1) = planet->m;
            display.block(0,3,J,4) = planet->y ;
            std::cout << display << std::endl;
            
            std::cout << planet->M << " scale " << planet->scale.transpose() << " Jcore " << planet->Jcore << std::endl ;

            //static_struct(planet , eos , acc);
            }
            
            static_struct(planet , eos , acc);
            
            std::cout << "Iron Core introduced.  Jcore " << planet->Jcore << std::endl;
            print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit);
            print.save_values(*num , planet , eos  , time , f_time , con);
            planet->nabtype = nabtype;

        }
        

        void Lcore_drop(struct P_P *planet , EoS *eos , double Lcore ,  double acc ){
            /** Function that removes imposed core luminosity so timestepping can start */
            
            double Lstep ;
            Linit = planet->scale(3);
          con.nmax = 2000;
		f_time = 1e-4;
		while (planet->L_core > Lcore){
		       
		Lstep = Linit*0.1;
		
                if (planet->L_core < Linit*0.2){    
		            Lstep = 0.2*planet->L_core;
                }
	            if(planet->L_core < 0.01*Linit){
		            
		            Lstep = planet->L_core-Lcore;
	            }else if(planet->L_core < - 2e19){
		            f_time = 0.001;
                    Lstep = planet->scale(3)*0.05;
	            }else{
		            //f_time = 0.00001;
		           // Lstep = planet->scale(3)*0.01;
	            }
	 if(Lstep >= planet->L_core - Lcore){
			        Lstep = planet->L_core - Lcore;
		      }	
		planet->L_core -= Lstep;
	        step_chooser(planet, f_time);
                planet->T0 = planet->scale(2) * planet->y.col(2);
                planet->P0 = planet->scale(0) * planet->y.col(0);
                if (planet->Jcore >0){
                    if( planet->phi(planet->Jcore) > 0.39){
                        planet->R_core0 = 0;
                    }else{
                        planet->R_core0 = planet->ym(planet->Jcore-1,1);
                    }
                }else{
                    planet->R_core0 = planet->R_core;
                }
                planet->R0 = planet->scale(1);
                con.converge(planet , eos , acc);
                std::cout << "Lcore " <<  planet->L_core << std::endl;
		        print.print_accuracies(planet , t , time , con );
            	if(con.maxinacc > acc){
		            std::cout << "\nHigh error!" <<std::endl;
			        if(con.maxinacc > 10*acc){
                        print.Not_Converged(planet , eos , con);
                        }else{
			            std::cout << "but accepted becaue not that high" << std::endl;
		            }
                }
     
		   time = time + planet->dt;
                
            }

        }
        void H_drop(struct P_P *planet , EoS *eos , double H ,  double acc ){
            /** Function that removes imposed internal heating so timestepping can start */

            double H0 = planet->H(50);
            //planet->H = planet->H - Hstep;
            f_time = 1e-3;
            if(planet->RadioType == 'E'){
                planet->H.tail(J-1-planet->Jcore) = radio.Earth(J-1-planet->Jcore,time/yrs_to_s);
                if(planet->Jcore > 0){
                    planet->H.head(planet->Jcore) = 0;
                }
            }else if (planet->RadioType == 'c'){
		        Eigen::ArrayXd Hstep = planet->H - H;		             
                if((planet->H>0).all() ){
                    planet->H = planet->H - Hstep;
                }
                if((planet->H < 0).any())
                {
                    planet->H = ArrayXd::Zero(J-1);
                } 
                radio.set_Constant(H);

            }
            step_chooser(planet, f_time);
            planet->T0 = planet->scale(2) * planet->y.col(2);
            planet->P0 = planet->scale(0) * planet->y.col(0);

            if (planet->Jcore >0){
                if( planet->phi(planet->Jcore) > 0.39){
                    planet->R_core0 = 0;
                }else{
                    planet->R_core0 = planet->ym(planet->Jcore-1,1);
                }
            }else{
                planet->R_core0 = planet->R_core;
            }   
            planet->R0 = planet->scale(1);

            con.converge(planet , eos , acc);
            if(con.maxinacc > acc){    
                print.save_values(1 , planet, eos , time , f_time , con);
                print.Not_Converged(planet , eos , con);
       	    }
            print.print_accuracies(planet , t , time , con );
            
        }
        void timestepper_Hdecay(struct P_P *planet , EoS *eos , int tmax , double yrs_max , double Lcore ,  double acc ,  double _f_time , int min_time_power , char *num  ,  std::string &_folder){
            /** Timestep assuming that the internal heating, H decays exponentially */

            t=0;
            time = 0;
            double time0;
            int time_power;
            bool save;
            
            ArrayXd orig = planet->scale , H0 = planet->H;
            planet->L_core = Lcore;
            
            con.nmax = 500;
	        f_time = _f_time;

	        print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit); 
            static_struct(planet , eos , acc);
            print.save_values(*num , planet , eos  , time , f_time , con);
            print.save_bulk(planet , eos , t , time );
	        f_time = _f_time;

            for(t=1; t<=tmax ; t++){
                save = false;
		        time0 = time;
		        if(t>500){
			        f_time = _f_time;
		        }
		        if(t>5000){
			        f_time = _f_time;
		        }
	            std::cout << "ftime " << f_time << std::endl;
                 
                step_chooser(planet, f_time);
                std::cout << "dt " << planet->dt << std::endl;
                time_power = std::ceil(log10( time0*1.000001/yrs_to_s )) ;
                time_power = std::max(time_power, 0);
                time = time0 + planet->dt;
                //std::cout << "time " << time/yrs_to_s << " time0 " <<  time0/yrs_to_s << " POWER " << time_power << " log " << log10( time/yrs_to_s )  << std::endl;
                if(log10( time/yrs_to_s ) > time_power && log10( time/yrs_to_s ) > min_time_power){
                    planet->dt = pow(10,time_power) * yrs_to_s - time0;
                    time = time0 + planet->dt;
                    save = true;
                    std::cout << "new time " << time/yrs_to_s << std::endl;
                }

                planet->H = H0*exp(-time/(1e4*yrs_to_s));
                std::cout << "H " << planet->H(0)/H0(0);
                planet->T0 = planet->scale(2) * planet->y.col(2);
                planet->P0 = planet->scale(0) * planet->y.col(0);
                if (planet->Jcore >0){
                    if( planet->phi(planet->Jcore) > 0.39){
                        planet->R_core0 = 0;
                    }else{
                        planet->R_core0 = planet->ym(planet->Jcore-1,1);
                    }
                }else{
                    planet->R_core0 = planet->R_core;
                } 
                        planet->R0 = planet->scale(1);

                con.converge(planet , eos , acc);

                if( (((planet->y.col(0)*planet->scale(0)-planet->P0)/planet->P0).abs()> 0.1).any() == 1 ){
           //     std::cout  << "P " << ((planet->y.block(0,0,J-1,1)*planet->scale(0)-planet->P0)/planet->P0).transpose() <<"\n";
                    
                   std::cout << "WOULD Reducing timestep because P drop too large" << std::endl;
               //     f_time = 0.5*f_time;
                 //   step_chooser(planet, f_time); 
                   // con.converge(planet , acc);
                }

                std::cout << t << " / " << tmax << " timestep (yrs) " << planet->dt/(3600.0*24*365) << std::endl;
                
                std::cout << "P / P_orig" << planet->scale(0)/orig(0) <<  " R/ R_orig" << planet->scale(1)/orig(1) << "T / T_orig" << planet->scale(2)/orig(2) <<  " L/ L_orig" << planet->scale(3)/orig(3) << std::endl;
                
                
                print.print_accuracies(planet , t , time , con );
                if( *num != '0'){
                     print.save_bulk(planet , eos , t , time);
                }
                if(con.maxinacc > acc){
                    std::cout << "\nHigh error!" <<std::endl;
                    if(con.maxinacc > 5*acc){
                        print.Not_Converged(planet , eos , con);
       	            }else{
                          std::cout << "but accepted becaue not that high" << std::endl;
                    }
             	}  
                if (save == true || t == tmax || con.maxinacc > acc){
                    print.save_values(*num , planet , eos  , time , f_time , con);
                }
                if (time/yrs_to_s >= yrs_max){
                    print.save_values(*num , planet , eos  , time , f_time , con);
                    break;
                }
            }

        }

        void static_mass_loss(struct P_P *planet , EoS *eos , int tmax , double yrs_max ,  double acc , int min_time_power , char *num  , char *log_or_lin , std::string &_folder){
            /** Timestep with no thermal evolution, only mass loss */

            std::cout << "running the static" << std::endl;
            t=0;
            time = 0;
            f_time = 1;
            double time0 , time_next , M0 , Mtot0 ;
            double time_power = min_time_power ;
            double min_lmix;
            ArrayXd orig = planet->scale , m0 , r0 , dm0;

            int j_fixed = int((1-0.05)*J) , J0;

            if(*log_or_lin == '1'){
                time_next = min_time_power*yrs_to_s;
            }else
            if(*log_or_lin == '2'){
                time_next = pow(10 , min_time_power)*yrs_to_s;
            }

            planet->M_total = planet->M;
            planet->M_edge = 0;
            planet->z_edge = 0;

            planet->M_from_day = 0;

            print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit); 
            //f_time = 1e-3 ; 
            print.save_values(*num , planet , eos  , time , f_time , con);
            print.save_bulk(planet , eos , t , time );
            print.save_composition(*num , planet , time , f_time);
            
            M_loss.read_table(planet->MassLossType , planet->Teq);
            calc_Mdot(planet);

            for(t=1; t<=tmax ; t++){
                J0 = J;
		        time0 = time;
                m0 = planet->m;
                M0 = planet->M;
                Mtot0 = planet->M_total;

                if (planet->Jcore >0){
                
                    //adjust_ML++;
                    min_lmix  = 0.05*planet->ym(planet->Jcore,0)/(planet->rho(planet->Jcore)*planet->g(planet->Jcore));
                    std::cout << "min mix " << min_lmix << std::endl;
                    planet->R_core0 = planet->scale(1)*planet->y(planet->Jcore-1,1) - min_lmix;
                    
                    std::cout << planet->R_core0 << " Rcore " << planet->R_core << std::endl;
                }else{
                    planet->R_core0 = planet->R_core;
                } 
                planet->R0 = planet->scale(1);
        
                j_fixed = int((1-0.05)*J) ;
            
                calc_Mdot(planet);

                if(planet->MassLossType == 'B' && planet->rho_fac < 0.7){
                    time = yrs_max*2.0;
                    print.save_bulk(planet , eos  , t , time );
                }
                std::cout << "mass-loss " << planet->M_dot <<  std::endl;
                if(planet->M_dot * planet->dt/planet->M >= 0.03  * planet->dm.tail(J - j_fixed).sum()){
                    planet->dt = planet->dm.tail(J - j_fixed).sum() * 0.03 * planet->M/planet->M_dot;
                    planet->dt = std::min(planet->dt , 2e8*yrs_to_s);
                    std::cout << "Timestep is mass-loss limited" << std::endl; 
                }
                planet->g_surf = G_Newt*planet->M/(planet->scale(1)*planet->scale(1));  
                
                time = time0 + planet->dt;
                std::cout << "first timestep " << planet->dt/yrs_to_s << std::endl;
                //std::cout << "time " << time/yrs_to_s << " time0 " <<  time0/yrs_to_s << " POWER " << time_power << " log " << log10( time/yrs_to_s )  << std::endl;
                std::cout << time/yrs_to_s << " next " << time_next/yrs_to_s << std::endl;
                
                if(time > time_next){
                    time = time0;
                    print.save_values(*num , planet , eos  , time , f_time , con);
                    time = time0 + planet->dt;
                }
                if (planet->MassLossType != '0'){
                    Merge_Regrid(planet , eos);
                }
                m0 = planet->m;
                r0 = planet->scale(1)*planet->y.col(1);
                dm0 = planet->dm;
                j_fixed = int((1-0.05)*J);

                std::cout << "M " << planet->M << " M total " << planet->M_total << std::endl;
                
                std::cout << "z edge " << planet->z_edge << " T " << planet->ym(J-1,2) << " Teq " << planet->Teq << std::endl;
                
                Lose_Mass( planet , planet->M_dot * planet->dt , j_fixed);
                
                planet->T0 = planet->scale(2) * planet->y.col(2);
                planet->P0 = planet->scale(0) * planet->y.col(0);
                planet->H = 0;
                planet->y.col(3) = 1;
                planet->L_core = planet->scale(3);
                planet->E_mass_loss = 0;
                planet->dt = 1e308;

                con.converge(planet , eos , acc);

                print.print_accuracies(planet , t , time , con );
                if( *num != '0'){
                     print.save_bulk(planet , eos , t , time );
                }
                if(con.maxinacc > acc){
                    std::cout << "H " << planet->E_mass_loss.transpose() << std::endl; 
                    print.save_values(*num , planet , eos  , time , f_time , con);
                    print.Not_Converged(planet , eos , con);
                    
       	        }
                if(time > time_next ){
                    std::cout << "onto next time output" << std::endl;

                    print.save_average_values(num , planet , eos , time_next , time , time0);
                    print.save_values(*num , planet , eos  , time , f_time , con);
                    
                    if(*log_or_lin== '1' ){
                        time_next = time_next + min_time_power * yrs_to_s;
                    }
                    else if (*log_or_lin == '2' /*&& time/yrs_to_s <= 1e8*/){
                        time_power = time_power+0.25 ;
                        time_next = std::ceil(pow(10,time_power)/10)*10 * yrs_to_s;
                    }
                    else if( time/yrs_to_s > 1e8){
                        time_next = time_next + 1e8*yrs_to_s;
                    }
                }else
                if (t%100 == 0|| t == tmax ){
                    std::cout << "100 steps " << std::endl;
                    print.save_values(*num , planet , eos  , time , f_time , con);
                } 
                if (time/yrs_to_s >= yrs_max || planet->M/M_Earth < 1e-7){
                    print.save_values(*num , planet , eos  , time , f_time , con);
                    break;
                }
            }
        }

        void timestepper(struct P_P *planet , EoS *eos , int tmax , double yrs_max , double Lcore ,  double acc ,  double _f_time , int min_time_power , char *num  , char *log_or_lin , std::string &_folder){
            /** Step forward in time for ascribed steps / time limit  */

            t=0;
            time = 0;
            double time0 , time_next , M0 , Mtot0 , dt0 , L0 , Lpred , Tspred , T_tol , L_tol , phi_c0 , phi_c1 ;
            double time_power = min_time_power ;
            double up_acc = 50*acc;
            double min_lmix;
            ArrayXd orig = planet->scale , m0 , r0 , dm0;

            int j_fixed = int((1-0.05)*J) , J0;
            int reduce , adjust_ML = 0;
            
            bool passed = false;

            if(*log_or_lin == '1'){
                time_next = min_time_power*yrs_to_s;
            }else
            if(*log_or_lin == '2'){
                time_next = pow(10 , min_time_power)*yrs_to_s;
            }
            
            planet->L_core = Lcore;
            planet->M_from_day = 0;

            con.nmax = 50;
	        f_time = _f_time;
            
	        print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit); 
            //f_time = 1e-3 ; 
            print.save_values(*num , planet , eos  , time , f_time , con);
            print.save_bulk(planet , eos , t , time );
            print.save_composition(*num , planet , time , f_time);
            
            for(t=1; t<=tmax ; t++){
                J0 = J;
		        time0 = time;
                dt0 = planet->dt;
                m0 = planet->m;
                M0 = planet->M;
                Mtot0 = planet->M_total;

                planet->solid_core_mass[0] = planet->solid_core_mass[1];
                planet->solid_core_mass[1] = planet->solid_core_mass[2];

                if( eos->EoStype == "Mix"){
                    planet->M_edge = 4*M_PI*planet->scale(1)*planet->scale(1)*(planet->P_BL+planet->Psurf[0])/planet->g(J-1);
                    planet->z_edge = planet->P_BL+planet->Psurf[0]/(planet->rho(J-1)*planet->g(J-1));
                }else{
                    planet->M_edge = 0;
                    planet->z_edge = 0;
                }
                if (planet->Jcore >0){
                std::cout << "phi inner " << planet->phi(planet->Jcore) << " phi outer " << planet->phi(J-1) << std::endl;
                std::cout << planet->solid_core_mass[0] << " core solids " << planet->solid_core_mass[1] <<  " " << planet->solid_core_mass[2] << std::endl;

                    if( planet->M > M_Earth*0.8){
                        planet->R_core0 = 0;
                    }else{
                        //adjust_ML++;
                        min_lmix  = 0.05*planet->ym(planet->Jcore,0)/(planet->rho(planet->Jcore)*planet->g(planet->Jcore));
                        std::cout << "min mix " << min_lmix << std::endl;
                        planet->R_core0 = planet->scale(1)*planet->y(planet->Jcore-1,1) - min_lmix;
                    }
                    std::cout << planet->R_core0 << std::endl;
                }else{
                    planet->R_core0 = planet->R_core;
                }
                planet->R0 = planet->scale(1);

		        if(t>500){
			        //f_time = _f_time;
		        }
		        if(t>5000){
			       // f_time = _f_time;
		        }
	            std::cout << "ftime " << f_time << std::endl;
                if(time/yrs_to_s<4e2 || eos->EoStype != "Mix"){
                    step_chooser(planet, f_time);
                }else{
                    passed = true;

                    if(planet->Tbound == 'm' && planet->scale(3) < 3.5e17 && planet->ym(J-1,2) > 1600){
                        L_tol = 0.1;
                    }else{
                        L_tol = 0.02;

                    }
                    if(planet->Tbound == 't' && planet->Tsurf[0] < 1690 && planet->Tsurf[0] > 1400){
                        T_tol = 0.0001 + std::max(0.0,(planet->Tsurf[0]-1655)/(1690-1655)*0.0019);
                        L_tol = 0.01;
                    }else if(planet->Tbound == 't' && planet->Tsurf[0] < 1810 && planet->Tsurf[0] > 1760){
                        T_tol = 0.002;
                        L_tol = 0.02;
                    }else if(planet->Tbound == 'b' && planet->BL == '1' && planet->ym(J-1,2) < 1645 && planet->ym(J-1,2) > 1600){
                        T_tol = 0.0001 + std::max(0.0,(planet->Tsurf[0]-1655)/(1690-1655)*0.0019);
                        L_tol = 0.01;
                    }else if(planet->Jcore > 0 && planet->phi(planet->Jcore) < 0.16 && planet->phi(planet->Jcore)> 0.0){
                        L_tol = 0.05;//pow((planet->phi(planet->Jcore)-0.1498)/0.01,2)*0.01 + 0.00001;
                        T_tol = 0.02;
                    }else if(planet->Jcore > 0 && planet->phi(planet->Jcore) > 0.35){
                        T_tol = 0.02;
                        L_tol = 0.01;
                    }else if(planet->Jcore > 0 && planet->phi(J-1) < 0.59 && planet->phi(J-1) > 0.39 && planet->Psurf[0]==1e8 && planet->Tbound == 'l'){
                        T_tol = 0.02;
                        L_tol = 0.1;
                    }else{
                        T_tol = 0.02;
                        L_tol = 0.3;
                    }
                    planet->dt = dt0/(abs(L0-planet->scale(3))/planet->scale(3))*L_tol;
                    Lpred = planet->scale(3) + planet->dt/dt0 * (planet->scale(3)-L0);
                    if(planet->Tbound == 't' || planet->Tbound == 'l'){

                        Tspred = con.Jacob.boundary.T_surf<double>(planet , Lpred);
                  
                        std::cout << Lpred/planet->scale(3) << " T pred " << Tspred << " Told " << planet->Tsurf[0] << std::endl;
                        planet->dt = std::min(1.0,T_tol/(abs(Tspred - planet->Tsurf[0])/planet->Tsurf[0]))*planet->dt;
                    }
                    std::cout << L_tol << " " << L_tol / (abs(L0-planet->scale(3))/planet->scale(3))  << " "  << T_tol/(abs(Tspred - planet->Tsurf[0])/planet->Tsurf[0]) << " adjusting " << planet->dt/dt0 << std::endl;


                }
                L0 = planet->scale(3);
                
                phi_c1 = phi_c0;
                phi_c0 = planet->phi(planet->Jcore);
                std::cout << phi_c1 << " phi_c " << phi_c0 << std::endl;

                std::cout << "dt " << planet->dt << std::endl;
                
                calc_Mdot(planet);
                if(planet->MassLossType == 'B' && planet->rho_fac < 0.7){
                    time = yrs_max*2.0;
                    print.save_bulk(planet , eos  , t , time );
                }
                if (planet->M_dot != 0 ){
                    j_fixed = int((1-0.05)*J) ;
                    std::cout << "mass-loss " << planet->M_dot <<  std::endl;
                    if(planet->M_dot * planet->dt/planet->M >= 0.03  * planet->dm.tail(J - j_fixed).sum()){
                        planet->dt = planet->dm.tail(J - j_fixed).sum() * 0.03 * planet->M/planet->M_dot;
                        std::cout << "Timestep is mass-loss limited" << std::endl; 
                        std::cout << "dt " << planet->dt << std::endl;
                    }
                    
                }
                if(planet->phi(planet->Jcore) >= 0.4 ){
                    planet->dt = std::min(6*yrs_to_s , planet->dt);
                }
               
                
                if(planet->dt > dt0*1.2/* && abs(planet->phi(planet->Jcore)-0.4) > 0.01*/){
                    std::cout << "not too much bigger" << std ::endl;
                    planet->dt = dt0*1.2;
                }
                if(planet->phi(J-1) < 0.3 && planet->phi(J-1) < 0.3  ){
                    planet->dt = std::max(50*yrs_to_s , planet->dt);
                }
        
                planet->g_surf = G_Newt*planet->M/(planet->scale(1)*planet->scale(1));
                
                time = time0 + planet->dt;
                std::cout << "first timestep " << planet->dt/yrs_to_s << std::endl;
                //std::cout << "time " << time/yrs_to_s << " time0 " <<  time0/yrs_to_s << " POWER " << time_power << " log " << log10( time/yrs_to_s )  << std::endl;
                std::cout << time/yrs_to_s << " next " << time_next/yrs_to_s << std::endl;
                
                if(time > time_next){
                    time = time0;
                    print.save_values(*num , planet , eos  , time , f_time , con);
                    time = time0 + planet->dt;
                }
             //   if(log10( time/yrs_to_s ) > time_power && log10( time/yrs_to_s ) > min_time_power){
               //     planet->dt = pow(10,time_power) * yrs_to_s - time0;
                 //   time = time0 + planet->dt;
                    
                //}
                
                
                std::cout << planet->m(j_fixed) << " m J " << planet->m(J-1) << std::endl;
                
                if (planet->MassLossType != '0'){
                    Merge_Regrid(planet , eos);
                }

                m0 = planet->m;
                r0 = planet->scale(1)*planet->y.col(1);
                dm0 = planet->dm;
                j_fixed = int((1-0.05)*J);

                std::cout << "M " << planet->M << " M total " << planet->M_total << std::endl;
                
                std::cout << "z edge " << planet->z_edge << " T " << planet->ym(J-1,2) << " Teq " << planet->Teq << std::endl;
                std::cout << "mdot " << planet-> M_dot << std::endl;
                Lose_Mass( planet , planet->M_dot * planet->dt , j_fixed);


                if(planet->M_dot != 0 ){
                    planet->E_mass_loss = Mass_Term(planet , m0 , r0 , dm0 , M0 , j_fixed , J);

                }else{
                    planet->M = planet->M_total - planet->M_edge;
                    planet->E_mass_loss = 0;
                }
                //std::cout << planet->y.col(3).transpose() << std::endl;

                planet->T0 = planet->scale(2) * planet->y.col(2);
                planet->P0 = planet->scale(0) * planet->y.col(0);
                
                Radioactivity_and_core(planet , time , dt0);

                if( planet->L_conv(planet->Jcore) == 0 && passed == true){
                    std::cout << "row now " << planet->y.row(0) << std::endl;
                    passed = true;
                }

                con.converge(planet , eos , acc);

                /*if(abs(planet->solid_core_mass[2]-planet->solid_core_mass[1])/planet->M_core> 0.01){
                    std::cout << "solidified too much " << planet->solid_core_mass[2]-planet->solid_core_mass[1] << std::endl;
                    con.maxinacc = up_acc*1000;
                    planet->dt = planet->dt * 0.01/(abs(planet->solid_core_mass[2]-planet->solid_core_mass[1])/planet->M_core);
                }*/
                if(con.maxinacc > up_acc){
                if( con.maxinacc > up_acc && phi_c0 < 0.06 && phi_c0 > 0 ){
                    
                    planet->dt = std::min(dt0*10 , abs(dt0 * phi_c0/(phi_c1-phi_c0)));

                    time = time0 + planet->dt;
                    planet->M_total = Mtot0;
                    planet->M = M0;
                    planet->m = m0;
                    Lose_Mass( planet , planet->M_dot * planet->dt , j_fixed);
                    Radioactivity_and_core(planet , time , dt0);

                    if(planet->M_dot != 0 ){                        
                        planet->E_mass_loss = Mass_Term(planet , m0 , r0 , dm0 , M0 , j_fixed , J);
                    }
                    std::cout << "Trying to skip the tricky low phi bit: dt " << planet->dt/yrs_to_s << " yrs " << std::endl;
                    planet->scale(2) = planet->scale(2)*0.9999;
                    con.converge(planet , eos , acc);

                }
                if( con.maxinacc > up_acc && abs(phi_c0-0.4) < 0.01 && phi_c0 > 0.4){
                    
                    planet->dt = 200*yrs_to_s;
                    time = time0 + planet->dt;
                    planet->M_total = Mtot0;
                    planet->M = M0;
                    planet->m = m0;
                    Lose_Mass( planet , planet->M_dot * planet->dt , j_fixed);
                    Radioactivity_and_core(planet , time , dt0);

                    if(planet->M_dot != 0 ){                        
                        planet->E_mass_loss = Mass_Term(planet , m0 , r0 , dm0 , M0 , j_fixed , J);
                    }
                    std::cout << "Trying to skip the tricky phi=0.4 bit" << std::endl;
                    planet->scale(2) = planet->scale(2)*0.999;
                    con.converge(planet , eos , acc);

                }
                reduce = 0;
                while(con.maxinacc > up_acc  && reduce < 4){
                   
                    reduce++;
                    planet->dt = planet->dt*0.99;
                    time = time0 + planet->dt;
                    
                    planet->M_total = Mtot0;
                    planet->M = M0;

                    planet->m = m0;
                    Lose_Mass( planet , planet->M_dot * planet->dt , j_fixed);
                    Radioactivity_and_core(planet , time , dt0);
                
                    if(planet->M_dot != 0 ){                        
                        planet->E_mass_loss = Mass_Term(planet , m0 , r0 , dm0 , M0 , j_fixed , J);
                    }
                    std::cout << "\nReduced timestep " << reduce << " times" << std::endl;

                    con.converge(planet , eos , acc);

                }
                if(con.maxinacc > up_acc ){
                    planet->dt = planet->dt/pow(0.99,reduce);
                }
                while(con.maxinacc > up_acc  && reduce < 6){
        
                    reduce++;
                    planet->dt = planet->dt*1.1;
                    time = time0 + planet->dt;
                    planet->M_total = Mtot0;
                    planet->M = M0;
                    planet->m = m0;
                    Lose_Mass( planet , planet->M_dot * planet->dt , j_fixed);
                    Radioactivity_and_core(planet , time , dt0);

                    if(planet->M_dot != 0 ){                        
                        planet->E_mass_loss = Mass_Term(planet , m0 , r0 , dm0 , M0 , j_fixed , J);
                    }
                    std::cout << "\nIncreased timestep " << reduce << " times" << std::endl;

                    con.converge(planet , eos , acc);

                }
                if(con.maxinacc > up_acc && planet->phi(planet->Jcore) < 0.41 && planet->phi(planet->Jcore) > 0 ){
                    std::cout << "trying a big one" << std::endl;
                    print.save_values(*num , planet , eos  , time , f_time , con);

                    if(planet->M < 0.6*M_E){
                        planet->dt = std::max(400*yrs_to_s , planet->dt);
                    }else{
                        planet->dt = std::max(20*yrs_to_s, planet->dt);
                    }
                    time = time0 + planet->dt;
                    planet->M_total = Mtot0;
                    planet->M = M0;
                    planet->m = m0;
                    planet->y.col(2) = planet->T0/planet->scale(2);
                    planet->y.col(0) = planet->P0/planet->scale(0);

                    Lose_Mass( planet , planet->M_dot * planet->dt , j_fixed);
                    Radioactivity_and_core(planet , time , dt0);

                    if(planet->M_dot != 0 ){                        
                        planet->E_mass_loss = Mass_Term(planet , m0 , r0 , dm0 , M0 , j_fixed , J);
                    }

                    con.converge(planet , eos , acc);
                    print.save_values(*num , planet , eos  , time , f_time , con);
                }
                if(con.maxinacc > up_acc && planet->phi(planet->Jcore) < 0.38 && planet->phi(planet->Jcore) > 0 ){
                    std::cout << "trying a really big one" << std::endl;
                    planet->dt = std::max(0.5*time , planet->dt);
                    time = time0 + planet->dt;
                    planet->M_total = Mtot0;
                    planet->M = M0;
                    planet->m = m0;
                    planet->y.col(2) = planet->T0/planet->scale(2);
                    planet->y.col(0) = planet->P0/planet->scale(0);

                    Lose_Mass( planet , planet->M_dot * planet->dt , j_fixed);
                    Radioactivity_and_core(planet , time , dt0);

                    if(planet->M_dot != 0 ){                        
                        planet->E_mass_loss = Mass_Term(planet , m0 , r0 , dm0 , M0 , j_fixed , J);
                    }

                    con.converge(planet , eos , acc);
                }
               
                }
                
                if( (((planet->y.col(0)*planet->scale(0)-planet->P0)/planet->P0).abs()> 0.1).any() == 1 ){
           //     std::cout  << "P " << ((planet->y.block(0,0,J-1,1)*planet->scale(0)-planet->P0)/planet->P0).transpose() <<"\n";
                    
                   std::cout << "WOULD Reducing timestep because P drop too large" << std::endl;
               //     f_time = 0.5*f_time;
                 //   step_chooser(planet, f_time); 
                   // con.converge(planet , acc);
                }
                if(planet->BL != '0' && planet->z_BL > 0.5*(planet->scale(1)*(1 - planet->y(J-3,1)))){
                    std::cout << planet->z_BL << " > " << 0.5*planet->scale(1)*(1 - planet->y(J-3,1)) << std::endl;
                    throw std::invalid_argument("Conductive layer too thick" );
                }
                //planet->L_core = core.Lcore(planet->scale(2)*planet->y(0,2) , planet->T0(0) , planet->dt);
                std::cout << t << " / " << tmax << " timestep (yrs) " << planet->dt/(3600.0*24*365) << " reduced " << reduce << std::endl;
                
                std::cout << "P / P_orig" << planet->scale(0)/orig(0) <<  " R/ R_orig" << planet->scale(1)/orig(1) << "T / T_orig" << planet->scale(2)/orig(2) <<  " L/ L_orig" << planet->scale(3)/orig(3) << std::endl;
                
                
                print.print_accuracies(planet , t , time , con );
                //std::cout << planet->y.col(3).transpose() << std::endl;
                if( *num != '0'){
                     print.save_bulk(planet , eos , t , time );
                     print.save_composition(*num , planet , time , f_time);
                }
                if(con.maxinacc > acc){
                    std::cout << "\nHigh error!" <<std::endl;
                    if(con.maxinacc > up_acc){
                        print.save_values(*num , planet , eos  , time , f_time , con);
                        print.Not_Converged(planet , eos , con);
                    
       	            }else{
                        std::cout << "but accepted becaue not that high" << std::endl;
                        print.save_values(*num , planet , eos  , time , f_time , con);

                    }
             	} 
                if(time > time_next ){
                    print.save_average_values(num , planet , eos , time_next , time , time0);
                    print.save_values(*num , planet , eos  , time , f_time , con);
                    
                    if(*log_or_lin== '1' ){
                        time_next = time_next + min_time_power * yrs_to_s;
                    }
                    else if (*log_or_lin == '2' /*&& time/yrs_to_s <= 1e8*/){
                        time_power = time_power+0.25 ;
                        time_next = std::ceil(pow(10,time_power)/10)*10 * yrs_to_s;
                    }
                    else if( time/yrs_to_s > 1e8){
                        time_next = time_next + 1e8*yrs_to_s;
                    }
                }else
                if (t%1000 == 0 || t == tmax || abs(planet->solid_core_mass[2]-planet->solid_core_mass[1])>0.0){
                    print.save_values(*num , planet , eos  , time , f_time , con);
                } 
                if (time/yrs_to_s >= yrs_max || planet->M/M_Earth < 1e-7){
                    print.save_values(*num , planet , eos  , time , f_time , con);
                    break;
                }
            }

        }
        void timestepper_full(struct P_P *planet , EoS *eos, int tmax, double yrs_max , double Lcore ,  double acc ,  double _f_time , int min_time_power , char *num  , char *log_or_lin , std::string &_folder){
            /** Remove imposed core luminosity and then Timesteps */

            t=0;
            time = 0;

            planet->M_dot = 0;
            planet->M_total = planet->M ;

            con.nmax = 2000;
            static_struct(planet , eos , acc);
	        print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit); 
	        print.save_values(*num , planet , eos  , time , f_time , con);
            
            Lcore_drop(planet , eos , Lcore , acc );
       
	        std::cout << "reset time" << std::endl;  
	    
            timestepper(planet , eos , tmax , yrs_max , planet->L_core , acc , _f_time , min_time_power, num , log_or_lin , _folder);    
        }
        
        void timestepper_fullH(struct P_P *planet , EoS *eos , int tmax , double yrs_max , double Hfinal ,  double acc ,  double _f_time , int min_time_power , char *num  , char *log_or_lin , std::string &_folder){
            /** Function that removes imposed internal heating and timesteps */
            /// Setup some quantities
            t=0;
            time = 0;
            con.nmax = 200;
            char nabtype = 'z';
            double M_dot = planet->M_dot;
            
            planet->M_dot = 0;
            planet->M_total = planet->M ;
            
            planet->g_surf = G_Newt*planet->M/(planet->scale(1)*planet->scale(1));
            

            // planet->M_edge = 4*M_PI*planet->scale(1)*planet->scale(1)*(planet->P_BL+planet->Psurf[0])/planet->g(J-1);
            // planet->z_edge = planet->P_BL+planet->Psurf[0]/(planet->rho(J-1)*planet->g(J-1));
            
            //planet->M = planet->M_total - planet->M_edge;
            //planet->H = planet->H*planet->M_total/planet->M;

            /// First Converge to a valid static structure, with an adiabatic T structure
            if (planet->nabtype == 'c' || planet->nabtype == 'g' || planet->nabtype == 'B'){
                nabtype = planet->nabtype;
                planet->nabtype = 'a';
            }
            static_struct(planet , eos , acc);
            std::cout << "Static adiabatic structure 1 found\n" << std::endl;

            planet->g_surf = G_Newt*planet->M/(planet->scale(1)*planet->scale(1));
            planet->M_edge = 4*M_PI*planet->scale(1)*planet->scale(1)*(planet->P_BL+planet->Psurf[0])/planet->g(J-1);

            planet->M = planet->M_total - planet->M_edge;  
              
            static_struct(planet , eos , acc);
            std::cout << "Static adiabatic structure 2 found\n" << std::endl;

            print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit); 
	        print.save_values(*num , planet , eos  , time , f_time , con);

            /// Now change to the desired T structure formulation (method of calculating nabla) and calculate another static structure     
            if (nabtype != 'z' ){                
                planet->nabtype = nabtype;

                std::cout << "Linit " << planet->scale(3)/(4*M_PI*planet->scale(1)*planet->scale(1)) << std::endl;
                int power = log10(planet->scale(3)/(4*M_PI*planet->scale(1)*planet->scale(1)) );
                Linit = planet->scale(3);
                planet->scale(3) = std::floor(planet->scale(3)/(4*M_PI*planet->scale(1)*planet->scale(1)) /pow(10,power))*pow(10,power) *(4*M_PI*planet->scale(1)*planet->scale(1)) ;
                planet->H = planet->H *planet->scale(3)/Linit;
                std::cout <<  "Linit " << planet->scale(3)/(4*M_PI*planet->scale(1)*planet->scale(1))  << std::endl;
                Linit = planet->scale(3)/(4*M_PI*planet->scale(1)*planet->scale(1));
                static_struct(planet , eos , acc);
                std::cout << "Static non-adiabatic structure found" << std::endl;
                std::cout << "nabtype = " << planet->nabtype << std::endl;

                print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit); 
    	        print.save_values(*num , planet , eos  , time , f_time , con);
            }
            /// Remove the internal heating that was used for the static stucture
            H_drop(planet , eos , Hfinal, acc );
            
	        std::cout << "reset time" << std::endl;  

            if (planet->MassLossType != '0'){
                   
                M_loss.read_table(planet->MassLossType , planet->Teq);
                calc_Mdot(planet);
            }

	        print.set_files(planet ,  *num  , f_time ,  _folder , eos->EoStype , Linit); 
            planet->g_surf = G_Newt*planet->M/(planet->scale(1)*planet->scale(1));

            timestepper(planet , eos , tmax , yrs_max , planet->L_core , acc , _f_time , min_time_power , num , log_or_lin , _folder);    
        }
        
        void calc_Mdot(struct P_P *planet ){
            /** Calculates the mass loss rate for the planet's conditions, under different assumptions encoded as `MassLossType` */
            if (planet->MassLossType != '0'){

            if(planet->MassLossType == 'c'){
                planet->M_dot = M_loss.Constant();
            }else if(planet->MassLossType == 'i'){
                    planet->M_dot = M_loss.Isothermal(planet->M_total/M_Earth)*M_Earth/(1e9*yrs_to_s);
            }
            else if(planet->MassLossType == 'P'){
                planet->M_dot = M_loss.Perez_PowerLaw(planet->M_total/M_Earth)*M_Earth/(1e9*yrs_to_s);
            }
            else if(planet->MassLossType == 'p'){
                planet->M_dot = M_loss.PBC_tabulated(planet->M_total/M_Earth)*M_Earth/(1e9*yrs_to_s);
            }else if(planet->MassLossType == 'B'){
                double R_fort;
                R_fort = Fortney_R(log10(planet->M/M_Earth) , 2.0/3);
                planet->rho_fac = std::min(pow(planet->scale(1)/R_Earth/R_fort , 3) , 1.2);
                                    
                if(planet->rho_fac > 0.7){
                    planet->M_dot = M_loss.Bicubic_tabulated(planet->M_total/M_Earth , (planet->scale(1)+planet->z_edge)/R_Earth , planet->rho_fac);
                }
            }else if(planet->MassLossType == 'K'){
                double period = 3600*24 ;
                planet->M_dot = M_loss.Kang_anal(planet->g(J-1) , (planet->scale(1)+planet->z_edge) , period);              
            }
            planet->M_dot_from_day = M_loss.Kang_terminator(planet->g(J-1) , (planet->scale(1)+planet->z_edge));
            std::cout << "mdot " << planet->M_dot << std::endl;
            }else{
                planet->M_dot = 0;
            }
        }

        void Lose_Mass( struct P_P *planet , double DeltaM , int j_fixed){
            /** Removes mass from the outer part of the planet */

            double M0 = planet->M;
            Eigen::ArrayXd mm0 = planet->mm;

            planet->M_total = planet->M_total - DeltaM;
            
            planet->M_from_day = planet->M_from_day + planet->dt * planet->M_dot_from_day;
        
            planet->M = planet->M_total - planet->M_edge;
            planet->m.tail(J-j_fixed)= 1.0 - (1.0 - M0/planet->M*planet->m(j_fixed))*(1.0-planet->m.tail(J-j_fixed ))/(1.0-planet->m(j_fixed)) ;
            planet->m.head(j_fixed) = planet->m.head(j_fixed)*M0/planet->M;
            
            JacobMaker::prepm(planet);
             
            for(int j= j_fixed ; j < J-1 ; j++){
                planet->y(j,1) = planet->y(j-1,1) + planet->M*planet->dm(j)/(4*M_PI*planet->y(j-1,1)*planet->y(j-1,1)*planet->rho(j)*pow(planet->scale(1),3));
            }
            planet->y(J-1,1) = planet->y(J-2,1); 
            planet->scale(1) = planet->scale(1) * planet->y(J-1,1);
            planet->y.col(1) = planet->y.col(1)/planet->y(J-1,1);
            //planet->scale(3) = 0.9*planet->scale(3);
            //planet->scale(2) = planet->scale(2)*planet->M/M0;
            //planet->scale(0) = planet->scale(0)*planet->M/M0;
            //planet->Psurf[0] = planet->rho(J-1)*G_Newt*planet->M/(planet->scale(1)*planet->scale(1))*planet->z_edge;

            std::cout << "M " << planet->M << " M total " << planet->M_total << std::endl;

       //std::cout << planet->m.transpose()*planet->M / 1e24<< std::endl;

        }

        void Shrink_Regrid(struct P_P *planet , int _J, int j_fixed){
            /** Regrid by shrinking all the cells (not used) */
            

            std::cout << "in here. Jnew " << _J << std::endl;
            std::cout << planet->m(j_fixed-1)<< " m " << planet->m(j_fixed) << std::endl;
            std::cout << planet->mm(j_fixed-1)<< " mm " << planet->mm(j_fixed) << std::endl;

            Eigen::ArrayXd temp_m = planet->m;
            Eigen::ArrayXXd temp_y = planet->y;
            std::vector<double> m_interp , m_interp2 , mm_interp , mm_interp2 , P_interp, r_interp , T_interp , L_interp;
            
            mm_interp.push_back(planet->mm(j_fixed-2));
            P_interp.push_back(planet->y(j_fixed-2,0));
            T_interp.push_back(planet->y(j_fixed-2,2));

            for( int j=0 ; j<= J-j_fixed ; j++){
                m_interp.push_back(planet->m(j_fixed-1+j));
                mm_interp.push_back(planet->mm(j_fixed-1+j));
                P_interp.push_back(planet->y(j_fixed-1+j,0));
                r_interp.push_back(planet->y(j_fixed-2+j,1));
                T_interp.push_back(planet->y(j_fixed-1+j,2));
                L_interp.push_back(planet->y(j_fixed-2+j,3));
                m_interp2 = m_interp;
                mm_interp2 = mm_interp;

            }

            planet->resize(_J);
            std::cout << planet->scale.transpose() << std::endl;

            J = _J;
            planet->m.head(j_fixed) = temp_m.head(j_fixed);
            planet->y.block(0,0,j_fixed,I) = temp_y.block(0,0,j_fixed,I);
            planet->m.tail(J-j_fixed+1) = temp_m.head(J-j_fixed+1)/temp_m(J-j_fixed)*(1-planet->m(j_fixed-1)) + planet->m(j_fixed-1) ;
            JacobMaker::prepm(planet);
            
            boost::math::interpolators::pchip<std::vector<double>> Pspline(std::move(mm_interp2) , std::move(P_interp));
            mm_interp2 = mm_interp;
            boost::math::interpolators::pchip<std::vector<double>> rspline(std::move(m_interp2) , std::move(r_interp));
            m_interp2 = m_interp;
            boost::math::interpolators::pchip<std::vector<double>> Tspline(std::move(mm_interp2) , std::move(T_interp));
            boost::math::interpolators::pchip<std::vector<double>> Lspline(std::move(m_interp2) , std::move(L_interp));
            std::cout << planet->m(j_fixed-1)<< " m " << planet->m(j_fixed) << std::endl;
            std::cout << planet->mm(j_fixed-1)<< " mm " << planet->mm(j_fixed) << std::endl;
            for( int j=0 ; j< J-j_fixed ; j++){
                planet->y(j_fixed+j,0) = Pspline(planet->mm(j_fixed+j));
                planet->y(j_fixed+j-1,1) = rspline(planet->m(j_fixed+j));
                planet->y(j_fixed+j,2) = Tspline(planet->mm(j_fixed+j));
                planet->y(j_fixed+j-1,3) = Lspline(planet->m(j_fixed+j));
            }
            planet->y(J-1,1) = planet->y(J-2,1);
            planet->y(J-1,3) = planet->y(J-2,3);

            std::cout << "J " << J << std::endl;
            //std::cout << planet->m.transpose()  << std::endl;

            std::cout << planet->scale.transpose()  << std::endl;

            con.resize(J);
            print.resize(J);
        }

        void Merge_Regrid(struct P_P *planet , struct  EoS *eos ){
            /** Regrid by merging and splitting cells (used) */

            Eigen::Vector4d DY , change ;
            std::vector<double> m_new , P_new , R_new , T_new , L_new , m_interp , y_interp ;
            ArrayXd m0;
            ArrayXXd y0;

            double f , f2 ;
            int stop = 0;
            std::cout << "Core0 " << planet->R_core0 <<  std::endl;
            //std::cout << planet->l_mix.segment(planet->Jcore-1,20).transpose() << std::endl;
            std::cout << "shrinking! J " << J <<  std::endl;
            std::cout <<"1st row" << planet->y.row(0) << std::endl;

            m_new.push_back(planet->m(planet->Jcore));
            
            change(0) = 0.001;
            if(planet->scale(3) > 1e14 && planet->Psurf[0] == 1e9){
                change(1) = 0.01;
                change(2) = 0.01;
            }else if(planet->scale(3) > 1e14 && planet->Psurf[0] == 1e8){
                change(1) = 0.01;
                change(2) = 0.01;
            }else{
                change(1) = 0.1;
                change(2) = 0.5;
            }
            change(3) = 0.001;
            
            int j=planet->Jcore ;
            if(planet->Jcore > 0){
                m0 = planet->m.head(planet->Jcore);
                y0 = planet->y.block(0,0,planet->Jcore,4);
            }
            while(j<J-2){
                
                DY(0) = planet->dy(j,2) / planet->y(j+1,2);
                
               
                DY(1) = (planet->Cp(j+1)-planet->Cp(j))/planet->Cp(j+1) ; // std::max( abs(planet->l_mix(j)-planet->l_mix(std::max(j-1,planet->Jcore))) , abs(planet->l_mix(j+1)-planet->l_mix(j)))/planet->l_mix(j) ;
                
                if( planet->dy(j,0) > 1e-6){
                    DY(2) = planet->dy(j,0)/planet->y(j,0) ;
                }else{
                    DY(2) = 0 ;
                }
                DY(3) = planet->dm(j)/(1-planet->m(planet->Jcore));
                
                DY = DY.array()/change.array();
                f = DY.squaredNorm();
                
                f2 = 8;
                //std::cout << j << " Jcore " <<  planet->Jcore << " " << DY.transpose() << " f "  << f << " lmixes " << planet->l_mix(j+2) << " " <<planet->l_mix(j+1) << " " << planet->l_mix(j) << " dm " << planet->dm(j) << std::endl;
                if(j == J-3){
                    std::cout << "ends " << DY.transpose() << std::endl;
                }
                if((f > DY.size()*10 && planet->dm(j) > 0.0001 ) ){
                    std::cout << "splitting  " << j << " / " << J << " " << planet->m(j) << " " << planet->Cp(j) << " " << planet->Cp(j+1) <<  std::endl;
                    std::cout << DY.transpose() << " " << f <<  std::endl;
                    
                    m_new.push_back(planet->m(j) + 4.0*M_PI*planet->y(j-1,1)*planet->y(j-1,1)*0.5*(planet->y(j,1)-planet->y(j-1,1))*planet->rho(j)*pow(planet->scale(1),3)/planet->M);
                    std::cout << m_new[m_new.size()-1] << " " << (planet->m(j+1) + planet->m(j))/2 << std::endl;
                    if(j < J-3){ 
                        m_interp = {planet->m(j),planet->m(j+1),planet->m(j+2), planet->m(j+3)};
                        y_interp = {planet->y(j-1,1), planet->y(j,1), planet->y(j+1,1), planet->y(j+2,1)};
                    }else{
                        m_interp = {planet->m(j-1),planet->m(j),planet->m(j+1), planet->m(j+2)};
                        y_interp = {planet->y(j-2,1), planet->y(j-1,1), planet->y(j,1), planet->y(j+1,1)};       
                    }
                    boost::math::interpolators::pchip<std::vector<double>> Rspline(std::move(m_interp) , std::move(y_interp));
                    R_new.push_back(Rspline(m_new[m_new.size()-1]));
                    R_new.push_back(planet->y(j,1));

                    if(j < J-3){
                        m_interp = {planet->m(j),planet->m(j+1),planet->m(j+2), planet->m(j+3)};
                        y_interp = {planet->y(j-1,3), planet->y(j,3), planet->y(j+1,3), planet->y(j+2,3)};
                    }else{
                        m_interp = {planet->m(j-1),planet->m(j),planet->m(j+1), planet->m(j+2)};
                        y_interp = {planet->y(j-2,3), planet->y(j-2,3), planet->y(j,3), planet->y(j+1,3)};       
                    }
                    boost::math::interpolators::pchip<std::vector<double>> Lspline(std::move(m_interp) , std::move(y_interp));
                    L_new.push_back(Lspline(m_new[m_new.size()-1]));
                    L_new.push_back(planet->y(j,3));
                    
                   if(j == 0){
                        m_interp = {planet->m(0),planet->mm(j+1),planet->mm(j+2), planet->mm(j+3)};
                        y_interp = {planet->y(j,2), planet->y(j+1,2), planet->y(j+2,2), planet->y(j+3,2)};
                    }else{
                        m_interp = {planet->mm(j-1),planet->mm(j),planet->mm(j+1), planet->mm(j+2)};
                        y_interp = {planet->y(j-1,2), planet->y(j,2), planet->y(j+1,2), planet->y(j+2,2)};       
                    }
                    boost::math::interpolators::pchip<std::vector<double>> Tspline(std::move(m_interp) , std::move(y_interp));
                    if(j==0){
                        T_new.push_back(planet->y(0,2));
                    }else{
                        T_new.push_back(Tspline((m_new[m_new.size()-1]+planet->m(j))/2));
                    }
                    T_new.push_back(Tspline((m_new[m_new.size()-1]+planet->m(j+1))/2));
                    
                    if(j == 0){
                        m_interp = {planet->m(0),planet->mm(j+1),planet->mm(j+2), planet->mm(j+3)};
                        y_interp = {planet->y(j,0), planet->y(j+1,0), planet->y(j+2,0), planet->y(j+3,0)};
                    }else{
                        m_interp = {planet->mm(j-1),planet->mm(j),planet->mm(j+1), planet->mm(j+2)};
                        y_interp = {planet->y(j-1,0), planet->y(j,0), planet->y(j+1,0), planet->y(j+2,0)};       
                    }
                    boost::math::interpolators::pchip<std::vector<double>> Pspline(std::move(m_interp) , std::move(y_interp));
                    if(j==0){
                        P_new.push_back(planet->y(0,0));
                    }else{
                        P_new.push_back(Pspline((m_new[m_new.size()-1]+planet->m(j))/2));
                    }
                    P_new.push_back(Pspline((m_new[m_new.size()-1]+planet->m(j+1))/2));
                    std::cout << T_new[T_new.size()-3] << " Ts "<< T_new[T_new.size()-2] << " Ts " << T_new[T_new.size()-1] << std::endl;
                    m_new.push_back(planet->m(j+1));
            
                    std::cout << m_new[m_new.size()-3] << " ms "<< m_new[m_new.size()-2] << " ms " << m_new[m_new.size()-1] << " " << planet->m(j+1) << std::endl;
                    j++;

                }else
                if( f < DY.size()/2.0 && stop == 0 && j > planet->Jcore + 60 ){
                    std::cout << "might merge " << j << std::endl;
                    std::cout << DY.transpose() << " f "  << f << std::endl;
                    
                    DY(0) = planet->dy(j+1,2) / planet->y(j+2,2);
                    DY(1) = (planet->Cp(j+2)-planet->Cp(j+1))/planet->Cp(j+2) ;
                    //DY(1) = planet->dmm(j+1) / planet->mm(j+2);
                    if( planet->dy(j+1,0) > 1e-6){
                        DY(2) = planet->dy(j+1,0)/planet->y(j+1,0) ;
                    }else{
                        DY(2) = 0 ;
                    }

                    if(j == J-3){
                        std::cout << "the end one" << std::endl;
                    }
                    DY(3) = planet->dm(j+1)/(1-planet->m(planet->Jcore));

                    
                    DY = DY.array()/change.array();
                    f2 = DY.squaredNorm();

                }
                if( f + f2 < DY.size() && stop == 0 && planet->dmm(j+1) < 0.005 && j > planet->Jcore+1){
                    std::cout << "merging " << j << std::endl;
                    std::cout << DY.transpose() << std::endl;

                    m_new.push_back(planet->m(j+2));
                    R_new.push_back(planet->y(j+1,1));
                    L_new.push_back(planet->y(j+1,3));

                    if(j < J-3){
                        m_interp = {planet->mm(j),planet->mm(j+1),planet->mm(j+2), planet->mm(j+3)};
                        y_interp = {planet->y(j,2), planet->y(j+1,2), planet->y(j+2,2), planet->y(j+3,2)};
                    }else{
                        m_interp = {planet->mm(j-1),planet->mm(j),planet->mm(j+1), planet->mm(j+2)};
                        y_interp = {planet->y(j-1,2), planet->y(j,2), planet->y(j+1,2), planet->y(j+2,2)};       
                    }
                    boost::math::interpolators::pchip<std::vector<double>> Tspline(std::move(m_interp) , std::move(y_interp));
                    T_new.push_back(Tspline((planet->m(j+2)+planet->m(j))/2));
                    if(j < J-3){
                        m_interp = {planet->mm(j),planet->mm(j+1),planet->mm(j+2), planet->mm(j+3)};
                        y_interp = {planet->y(j,0), planet->y(j+1,0), planet->y(j+2,0), planet->y(j+3,0)};
                    }else{
                        m_interp = {planet->mm(j-1),planet->mm(j),planet->mm(j+1), planet->mm(j+2)};
                        y_interp = {planet->y(j-1,0), planet->y(j,0), planet->y(j+1,0), planet->y(j+2,0)};       
                    }
                    boost::math::interpolators::pchip<std::vector<double>> Pspline(std::move(m_interp) , std::move(y_interp));
                    P_new.push_back(Pspline((planet->m(j+2)+planet->m(j))/2));
                    j = j+2;
                    stop = 0;

                }else{
                    m_new.push_back(planet->m(j+1));
                    P_new.push_back(planet->y(j,0));
                    R_new.push_back(planet->y(j,1));
                    T_new.push_back(planet->y(j,2));
                    L_new.push_back(planet->y(j,3));
                    j++;
                }
            }
            std::cout << "run through " << j <<  std::endl;
            if( j < J-1){
                m_new.push_back(planet->m(J-1));
                P_new.push_back(planet->y(J-2,0));
                R_new.push_back(planet->y(J-2,1));
                T_new.push_back(planet->y(J-2,2));
                L_new.push_back(planet->y(J-2,3));
            }
            R_new.push_back(planet->y(J-1,1));
            L_new.push_back(planet->y(J-1,3));
            P_new.push_back(planet->y(J-1,0));
            T_new.push_back(planet->y(J-1,2));
            std::cout << P_new.size() << " " << P_new.size() << " " << R_new.size() << " " << T_new.size() << " " << L_new.size() << " size " << J << std::endl;
            if( m_new.size() != J - planet->Jcore || stop == 1){
                J = m_new.size()+planet->Jcore; 
                planet->resize(J);
                planet->m.head(planet->Jcore) = m0;
                planet->m.tail(J-planet->Jcore) = Eigen::Map< Eigen::ArrayXd>(m_new.data(),J-planet->Jcore);
                JacobMaker::prepm(planet);
               
                planet->y.block(0,0,planet->Jcore,4) = y0;
                planet->y.block(planet->Jcore,0,J-planet->Jcore,1) = Eigen::Map< Eigen::ArrayXd>(P_new.data(),J-planet->Jcore);
                planet->y.block(planet->Jcore,1,J-planet->Jcore,1)  = Eigen::Map< Eigen::ArrayXd>(R_new.data(),J-planet->Jcore);
                planet->y.block(planet->Jcore,2,J-planet->Jcore,1)  = Eigen::Map< Eigen::ArrayXd>(T_new.data(),J-planet->Jcore);
                planet->y.block(planet->Jcore,3,J-planet->Jcore,1)  = Eigen::Map< Eigen::ArrayXd>(L_new.data(),J-planet->Jcore);
                
                planet->ym.col(0)= planet->y.col(0)*planet->scale(0);
                planet->ym.col(2)= planet->y.col(2)*planet->scale(2);
                planet->P0 = planet->ym.col(0);
                planet->T0 = planet->ym.col(2);
                
                eos->compute_values<ArrayXd>(planet);

                planet->k_cond = Conduct<ArrayXd>::k(planet->rho , planet->ym.col(2), planet->conds[0],planet->conds[1],planet->conds[2]);
                con.resize(J);
                print.resize(J);
            }
            
            std::cout << "Resizing done" << std::endl;
            std::cout << planet->y.row(0) << std::endl;
            //std::cout << planet->m.segment(planet->Jcore-2, 30).transpose() << std::endl;
        }

        void static_struct(struct P_P *planet , EoS *eos ,  double acc){
            /** Compute the structure assuming no time evolution (t-dependent terms -> 0) */

            planet->dt = 1e308;
		    f_time = 0;
            planet->T0 = planet->scale(2) * planet->y.col(2);
            planet->P0 = planet->scale(0) * planet->y.col(0);
            if (planet->Jcore >0){
                
                if( planet->phi(planet->Jcore) > 0.39){
                    planet->R_core0 = 0;
                }else{
                    planet->R_core0 = planet->ym(planet->Jcore-1,1);
                }
                
            }else{
                planet->R_core0 = planet->R_core;
            } 
                        planet->R0 = planet->scale(1);

            con.converge(planet , eos , acc);
            print.print_accuracies(planet , t , time , con );
            

            if(con.maxinacc > acc){
                print.save_values('1' , planet, eos , time , f_time , con);
                print.Not_Converged(planet , eos , con);
            }
        }
        
        void step_chooser(struct P_P *planet , double f ){
            double Eth = planet->M*((planet->del_rho.head(J-1)*planet->scale(0)*planet->y.block(0,0,J-1,1)-planet->Cp.head(J-1)*planet->scale(2)*planet->y.block(0,2,J-1,1))*planet->dm).sum();
            planet->dt = abs(f*Eth /(planet->scale(3)));
            //if(  r_CRB/planet->scale(1) > 0.75){
              // planet->dt = (0.01 + pow(1-r_CRB/planet->scale(1),2)*2.2) * yrs_to_s ;
              // std::cout << "Close to regime change - shrink timestep" << std::endl;
            //}-4*M_PI*planet->scale(1)*planet->scale(1)*Sig_Stefan*pow(planet->Teq,4)
        }

    stepper(int _I, int _J):
        I(_I) , J(_J) , con(_I,_J) , print(_I,_J)
        {};
    
        Converger con;

    private:
        Radioactivity radio;
        Core_Latent core_latent;
        MassLoss<double> M_loss;
        int I , J , t;
        double time , f_time , Linit;
       
        void Radioactivity_and_core(struct P_P *planet , double &time , double &dt0){
            /** Computes radioactive heating, and heat contribution from an unresolved core (if it exists) */
            planet->solid_core_mass[2] = 0.0;
            if(planet->RadioType == 'c'){
                planet->H = radio.Constant(J-1);
            }else if(planet->RadioType == 'E'){
                planet->H.tail(J-1-planet->Jcore) = radio.Earth(J-1-planet->Jcore,time/yrs_to_s);
                planet->H.tail(J-1-planet->Jcore) = planet->H.tail(J-1-planet->Jcore) ;//+ radio.Al26_Earth(J-1-planet->Jcore,time/yrs_to_s);

                if(planet->Jcore > 0){
                    /// Do not generally use the latent heat from resolved core as it can cause issues
                    planet->H.head(planet->Jcore) = radio.Iron60_Earth(planet->Jcore,time/yrs_to_s); // + core_latent.Heat(planet->solid_core_mass , planet->M*planet->m(planet->Jcore) , planet->Jcore ,  dt0);
                }
            }else if (planet->RadioType == '0'){
                planet->H = 0.0;
            }
            if(planet->R_core > 0 && planet->Jcore == 0){
                Iron_core<double> core(planet->M_core);
                planet->scale(2) = planet->scale(2) + planet->L_core/core.dLcore_dT(planet->dt);
                    
            }
        }
        
       
};
