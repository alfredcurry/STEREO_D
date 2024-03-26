#include <iostream>
#include <Eigen/Dense>
#include "StructStruct.hpp"
#include "PhysicalConsts.h"
#include "converge.hpp"

using namespace Eigen;

class Saver{ /// Class for saving output values, and printing key values
    public:
        std::string bulkfile , massfile;
        
        void print_accuracies(struct P_P *planet , int &t , double &time , Converger &con ){
            //print out some key values
            std::cout << "M " << planet->M/M_Earth<< " M total " << planet->M_total/M_Earth << std::endl;
            if( planet->dt <1e308 || planet->M_dot > 0 ){
                std::cout << "t: " << t << " time: " << time/yrs_to_s << " yrs ";
            }
            if( planet->M_core != 0){
                std::cout << " M " << planet->M << " M_core " << planet->M_core << " R_core " << planet->R_core<< " L_core " << planet->L_core << std::endl;
            }
            std::cout << planet->scale.transpose() << " Psurf " << planet->Psurf[0] << " z edge " << planet->z_edge << " P_bl " << planet->ym(J-1,0) << " Tsurf " << planet->Tsurf[0] <<" T_bl " << planet->ym(J-1,2) ;
            
            if(planet->Tbound == 'i'){
                std::cout << " Teq " << planet->Teq ;
            }
            std::cout << std::endl;
            
            std::cout << "No. iterations: " << con.n << " of a maximum of: "<< con.nmax << std::endl;
             
            std::cout << "Max innacuracy: " << con.maxinacc << ", Mean Accuracy: " << con.meaninacc << "\n" << std::endl; 

        }

        void save_average_values(char *num , struct P_P *planet , EoS *eos , double &save_time , double &time , double &time0){
            /// Save values for a particular time, that the code does not exactly calculate for, by averaging the values from the step before and the step after
            if(Allthing.col(0).size() != J){
                // The model mass grid needs to be the same size for both timesteps 
                std::cout << "DIFFERENT SIZES SO HAVEN'T AVERAGED" << std::endl;
            }else{
            
            char timestring[100];
            int a;
            double avg_factor = (save_time - time0)/(time - time0);
          

            std::ofstream fout;
            a = snprintf(timestring, 100 , "t=%.5g.txt" , save_time/yrs_to_s);
            if( eos->EoStype == "Mix"){

                Allthing.col(I+4) = (1-avg_factor)*Allthing.col(I+4) + avg_factor*planet->phi;
                Allthing.block(planet->Jcore,I+5,J-planet->Jcore,1) = (1-avg_factor)*Allthing.block(planet->Jcore,I+5,J-planet->Jcore,1) + avg_factor*eos->melt.Tsol;
                Allthing.block(planet->Jcore,I+5,J-planet->Jcore,1) = (1-avg_factor)*Allthing.block(planet->Jcore,I+5,J-planet->Jcore,1) + avg_factor*eos->melt.Tliq;
                Allthing.block(0,I+7,J,4) = 0*Allthing.block(0,I+7,J,4)  ; 


                Fls.col(3) = (1-avg_factor)*Fls.col(3) + avg_factor*planet->L_grav;
                Fls.col(4) = (1-avg_factor)*Fls.col(4) + avg_factor*planet->L_mix;

                others.col(0) = (1-avg_factor)*others.col(0) + avg_factor*planet->Cp ;
                others.col(1) = (1-avg_factor)*others.col(1) + avg_factor*planet->del_rho*planet->rho;
                others.col(2) = (1-avg_factor)*others.col(2)  + avg_factor*planet->nu;
                others.col(3) = (1-avg_factor)*others.col(3) + avg_factor*planet->l_mix;
                others.col(4) = (1-avg_factor)*others.col(4) + avg_factor*planet->DeltaS;
                others.col(5) = (1-avg_factor)*others.col(5) + avg_factor*planet->DeltaV;
                others.col(6) = (1-avg_factor)*others.col(6) + avg_factor*planet->u_conv;
                others.col(7) = (1-avg_factor)*others.col(7) + avg_factor*planet->g;
                
                

                fout.open(folder+"/AVG"+"/other"+num+basicfile+timestring);
                        fout << others;
                fout.close();
              
            }
            Allthing.col(0) =  (1-avg_factor)*Allthing.col(0) + avg_factor*planet->m;
            Allthing.block(0,1,J,I) =  (1-avg_factor)*Allthing.block(0,1,J,I) + avg_factor*planet->y;
            Allthing.col(I+1) =  (1-avg_factor)*Allthing.col(I+1) + avg_factor*planet->rho;
            Allthing.col(I+2) =  (1-avg_factor)*Allthing.col(I+2) + avg_factor*planet->nabla;
            Allthing.col(I+3) =  (1-avg_factor)*Allthing.col(I+3) + avg_factor*planet->nabNonAd;

            Fls.col(0) = (1-avg_factor)*Fls.col(0) + avg_factor*planet->ym.col(3);
            Fls.col(1) = (1-avg_factor)*Fls.col(1) + avg_factor*planet->L_cond;
            Fls.col(2) = (1-avg_factor)*Fls.col(2) + avg_factor*planet->L_conv;

            if (*num == '0'){
                    std::cout << "no file saved" <<std::endl;
            }else{
                   
                    fout.open(folder+"/AVG"+"/results"+num+basicfile+timestring);
                    std::cout << folder+"/AVG"+"/results"+num+basicfile+timestring << std::endl;
                    //fout << "m/M $p/P_C$ r/R rho G1 G2 \n";
                        fout << Allthing ;
                    fout.close();
                    fout.open(folder+"/AVG"+"/scale"+num+basicfile+timestring);
                        fout << planet->M <<"\n"<< planet->scale ;
                        if(planet->M_core != 0){
                            fout <<"\n"<< planet->M_core << "\n" << planet->R_core0;
                        }
                    fout.close();
                    fout.open(folder+"/AVG"+"/Fluxes"+num+basicfile+timestring);
                        fout << Fls;
                    fout.close();
                
                    std::cout << "Average files saved\n" <<std::endl;
            }
        }}
        void save_composition(char num , struct P_P *planet , double &time , double &f_time){
            /// Save the internal values of key parameters that may be important for calculating compostional evolution afterwards
            char timestring[100];
            std::ofstream fout;
            int a ;
            if(f_time == 0){
                a = snprintf(timestring, 100 , ".txt" );
            }else if( planet->M_dot == 0){
                a = snprintf(timestring, 100 , "t=%.4g.txt" , time/yrs_to_s);
            }else{
                a = snprintf(timestring, 100 , "Mf=%.4gt=%.4g.txt" , planet->M_total/M_Earth/Minit , time/yrs_to_s);
            }
        
                composition_array = Eigen::ArrayXXd::Zero(J-planet->Jcore-1,7);
                composition_array.col(0) = planet->mm.tail(J-planet->Jcore-1)*planet->M;
                composition_array.col(1) = planet->ym.col(1).segment(planet->Jcore, J-planet->Jcore-1);
                composition_array.col(2) = planet->ym.col(0).segment(planet->Jcore, J-planet->Jcore-1);
                composition_array.col(3) = planet->ym.col(2).segment(planet->Jcore, J-planet->Jcore-1);
                composition_array.col(4) = planet->g.segment(planet->Jcore, J-planet->Jcore-1);
                composition_array.col(5) = planet->u_conv.segment(planet->Jcore, J-planet->Jcore-1)*planet->l_mix.segment(planet->Jcore, J-planet->Jcore-1);
                composition_array.col(6) = planet->phi.segment(planet->Jcore, J-planet->Jcore-1);

            if (time/yrs_to_s <= 1e-3){
                fout.open(folder+"/time"+num+basicfile+".txt" );
                    fout << time/yrs_to_s << std::endl;
                fout.close();
            }else{
                fout.precision(8);
                fout.open(folder+"/time"+num+basicfile+".txt", std::ios::app);
                    fout  << time/yrs_to_s << std::endl;
                fout.close();
            }
                
            if (num == '0'){
                    std::cout << "no file saved" <<std::endl;
            }else{
                   
                    fout.open(folder+"/comp"+num+basicfile+timestring);
                    std::cout << folder+"/comp"+num+basicfile+timestring << std::endl;
                    //fout << "m/M $p/P_C$ r/R rho G1 G2 \n";
                        fout << composition_array ;
                    fout.close();

                    std::cout << "comp file saved\n" <<std::endl;
            }
        }
        void save_values(char num , struct P_P *planet , EoS *eos , double &time , double &f_time , Converger &con ){
            /// Save the internal structure
            char timestring[100];

            int a ;
            if(f_time == 0){
                a = snprintf(timestring, 100 , ".txt" );
            }else if( planet->M_dot == 0){
                a = snprintf(timestring, 100 , "t=%.5g.txt" , time/yrs_to_s);
            }else{
                a = snprintf(timestring, 100 , "Mf=%.5gt=%.5g.txt" , planet->M_total/M_Earth/Minit , time/yrs_to_s);

            }

            std::ofstream fout;
            /// 'Allthing' contains most of the key parameters, plus the errors on P,r,T,L
            /// Fls contains the contribtions to the heat fluxes
            /// 'others' contains some more quantities

            if( eos->EoStype == "Mix"){
                Allthing = Eigen::ArrayXXd::Zero(J,I + I + 7);

                Fls = Eigen::ArrayXXd::Zero(J,5);

                Allthing.col(I+4) = planet->phi;
                Allthing.block(planet->Jcore,I+5,J-planet->Jcore,1) = eos->melt.Tsol;
                Allthing.block(planet->Jcore,I+6,J-planet->Jcore,1) = eos->melt.Tliq;
                Allthing.block(0,I+7,J,4) = con.ErrorArray ; 

                Fls.col(3) = planet->L_grav;
                Fls.col(4) = planet->L_mix;

                others = Eigen::ArrayXXd::Zero(J,8);
                others.col(0) = planet->Cp;
                others.col(1) = planet->del_rho*planet->rho;
                others.col(2) = planet->nu;
                others.col(3) = planet->l_mix;
                others.col(4) = planet->DeltaS;
                others.col(5) = planet->DeltaV;
                others.col(6) = planet->u_conv;
                others.col(7) = planet->g;


                fout.open(folder+"/other"+num+basicfile+timestring);
                        fout << others;
                fout.close();
            }else{

                Allthing = Eigen::ArrayXXd::Zero(J,I + I + 3);
                Fls = Eigen::ArrayXXd::Zero(J,3);
                
                Allthing.block(0,I+3, J,I) = con.ErrorArray ; 
            }

            Allthing.col(0) = planet->m;
            Allthing.block(0,1,J,I) = planet->y;
            Allthing.block(0, I+1,J,1) = planet->rho;
            Allthing.block(0,I+2,J,1) = planet->nabla;
            Allthing.block(0,I+3,J,1) = planet->nabNonAd;

            Fls.col(0) = planet->ym.col(3);
            Fls.col(1) = planet->L_cond;
            Fls.col(2) = planet->L_conv;

            /// 'debug_quantities' is for any extra quantities that might be useful for debugging
            debug_quantities = Eigen::ArrayXXd::Zero(J,4);
                
            debug_quantities.col(0) = - planet->Cp * (planet->scale(2) * planet->y.block(0,2,J,1)-planet->T0)/planet->dt ;      
            debug_quantities.col(1).head(J-1) = planet->H;
            debug_quantities.col(2) = planet->del_rho * (planet->scale(0)*planet->y.block(0,0,J,1)-planet->P0)/planet->dt ;
            debug_quantities.col(3).head(J-1) = planet->E_mass_loss;

            if (num == '0'){
                    std::cout << "no file saved" <<std::endl;
            }else{
                   
                    fout.open(folder+"/results"+num+basicfile+timestring);
                    std::cout << folder+"/results"+num+basicfile+timestring << std::endl;
                    //fout << "m/M $p/P_C$ r/R rho G1 G2 \n";
                        fout << Allthing ;
                    fout.close();
                    fout.open(folder+"/scale"+num+basicfile+timestring);
                        fout << planet->M <<"\n"<< planet->scale ;
                        if(planet->M_core != 0){
                            fout <<"\n"<< planet->M_core << "\n" << planet->R_core ;
                        }
                    fout.close();
                    fout.open(folder+"/Fluxes"+num+basicfile+timestring);
                        fout << Fls;
                    fout.close();
                    fout.open(folder+"/debug"+num+basicfile+timestring);
                        fout << debug_quantities;
                    fout.close();
                    
                   
                    std::cout << "files saved\n" <<std::endl;
            }
        }
        void calc_E(struct P_P *planet ){
            /// Function for calculating the internal energy. Used for the Ideal Gas structures
		    Etherm = -planet->M*((-planet->Cp.head(J-1)*planet->ym.block(0,2,J-1,1) + planet->del_rho.head(J-1)*planet->ym.block(0,0,J-1,1))*planet->dm).sum();
                double Etherm2 = planet->M * (boltz_R/planet->mu_m*planet->ym.block(0,2,J-1,1)/planet->n * planet->dm).sum();
                Eg = - G_Newt*planet->M*planet->M * (planet->mm.head(J-1)*planet->dm/planet->ym.block(0,1,J-1,1)).sum();
                double Etot2 = (1-planet->n/3)*Eg + planet->n * 4*M_PI/3 * planet->Psurf[0] * pow(planet->scale(1),3);
                
                Etot_calc = Eg + Etherm;
	    }
        void save_bulk( struct P_P *planet , EoS *eos , int &t , double &time ){
            /// Save various bulk quantities, including the position of melt fronts
            std::ofstream fout;    

            if (t == 0){
                fout.open(massfile );
            }else{
                fout.open(massfile, std::ios::app);
                fout << std::endl;
            }
            fout.precision(8);
            fout << time/yrs_to_s << " " << planet->M_total << " " << planet->M_dot << " " << planet->M_dot_from_day << " " << planet->M_from_day << " " << planet->scale(1)+planet->z_edge;
            fout.close();   


            if (t == 0){
                fout.open(bulkfile );
            }
            else{
                fout.open(bulkfile, std::ios::app);
                fout << std::endl;
            }
            fout.precision(8);
            fout << time/yrs_to_s << " " << planet->scale.transpose() << " " << planet->Psurf[0] << " " << planet->Tsurf[0] << " " ;
                
            r_RCB = planet->R_core;
            r_CRB = planet->R_core;
            double m_RCB= 0, m_CRB = 0, rho_RCB =0 , Temp_rat , r_B ;
            int J = planet->m.size() , j;
            calc_E(planet);

            if( eos->EoStype == "Mix"){
                r_solup = planet->ym(J-1,1);
                P_solup = planet->ym(J-1,0);
                m_solup = planet->m(J-1);
                r_001up = planet->ym(J-1,1);

                if(planet->phi(J-1) > 0){
                    for(j=J-1 ; j>=planet->Jcore ; j--){
                        if(planet->phi(j) == 0){
                            r_soldown = planet->ym(j,1);
                            P_soldown = planet->ym(j,0);
                            m_soldown = planet->m(j);
                        
                            break;
                        }
                    }
                    if(j<=planet->Jcore){
                        r_soldown = planet->ym(planet->Jcore,1);
                        P_soldown = planet->ym(planet->Jcore,0);
                        m_soldown = planet->m(planet->Jcore);
                    }
                }else{
                    for(j=J-1 ; j>=planet->Jcore ; j--){
                        if(planet->phi(j) > 0){
                            r_solup = planet->ym(j,1);
                            P_solup = planet->ym(j,0);
                            m_solup = planet->m(j);

                            break;
                        }
                    }
                    if(j<=planet->Jcore){
                        r_solup = planet->ym(planet->Jcore,1);
                        P_solup = planet->ym(planet->Jcore,0);
                        m_solup = planet->m(planet->Jcore);
                        r_soldown = planet->ym(J-1,1);
                        P_soldown = planet->ym(J-1,0);
                        m_soldown = planet->m(J-1);
                    }else{
                        for(j=j ; j>=planet->Jcore ; j--){
                            if(planet->phi(j) == 0){
                                r_soldown = planet->ym(j,1);
                                P_soldown = planet->ym(j,0);
                                m_soldown = planet->m(j);
                            
                                break;
                            }
                        }
                        if(j<=planet->Jcore){
                            r_soldown = planet->ym(planet->Jcore,1);
                            P_soldown = planet->ym(planet->Jcore,0);
                            m_soldown = planet->m(planet->Jcore);
                        }
                    }
                    
                }
                for(j=planet->Jcore ; j<J ; j++){
                    if(planet->phi(j) == 1){
                        r_liq = planet->ym(j,1);
                        P_liq = planet->ym(j,0);
                        m_liq = planet->m(j);

                        break;
                    }
                }
                for(j=planet->Jcore ; j<J ; j++){
                    if(planet->phi(j) >= 0.4){
                        r_crit = planet->ym(j,1);
                        P_crit = planet->ym(j,0);
                        m_crit = planet->m(j);


                        break;
                    }
                }
                if(r_crit > planet->ym(J-1,1) || j >= J-1){
                    r_crit = planet->ym(J-1,1);
                    P_crit = planet->ym(J-1,0);
                    m_crit = planet->m(J-1);
                }
                double essentially_solid = 0.05;
                if(planet->phi(J-1) > essentially_solid){
                    for(j=J-1 ; j>=planet->Jcore ; j--){
                        if(planet->phi(j) < essentially_solid){
                            r_001down = planet->ym(j,1);

                            break;
                        }
                    }
                    if(j<=planet->Jcore){
                        r_001down = planet->ym(planet->Jcore,1);
                    }
                }else{
                    for(j=J-1 ; j>=planet->Jcore ; j--){
                        if(planet->phi(j) > essentially_solid){
                            r_001up = planet->ym(j,1);
                            break;
                        }
                    }
                    if(j<=planet->Jcore){
                        r_001up = planet->ym(planet->Jcore,1);
                        r_001down = planet->ym(J-1,1);

                    }else{
                        for(j=j ; j>=planet->Jcore ; j--){
                            if(planet->phi(j) < essentially_solid){
                                r_001down = planet->ym(j,1);
                                
                                break;
                            }
                        }
                        if(j<=planet->Jcore){
                            r_001down = planet->ym(planet->Jcore,1);
                            }
                    }
                    
                }
                
                if(r_liq > planet->ym(J-1,1)){
                    r_liq = planet->ym(J-1,1);
                    P_liq = planet->ym(J-1,0);
                    m_liq = planet->m(J-1);

                }
               

                for(j=J-1 ; j>=planet->Jcore ; j--){
                    if (planet->L_conv(j) > 0){
                        r_RCB = planet->ym(j,1);
                        break;    
                    }
                                    
                }
                if(j == 0){
                    r_RCB = 0;
                }
                for(j=J-1 ; j>=0 ; j--){
                    if (planet->ym(j,0) > 1e9){
                        T1Gpa = planet->ym(j,2);
                        break;    
                    }
                                    
                }
            }else{                

                r_B = G_Newt * planet->M_core/(boltz_R*planet->Teq/planet->mu_m);
                for(j=0 ; j<J ; j++){
                    if (planet->nabla(j) >= planet->nabAd(j)){
                        m_CRB = planet->m(j)  - planet->M_core / planet->M;
			            r_CRB = planet->ym(j,1);
                        break;
                    }
                }
                
                for(j=J-1 ; j>=0 ; j--){
                    if (planet->nabla(j) >= planet->nabAd(j)){
                        m_RCB = planet->m(j) - planet->M_core / planet->M ;
			            r_RCB = planet->ym(j,1);
                        rho_RCB = planet->rho(j);
                        Temp_rat = planet->y(j,2)/planet->y(J-1,2);
                        break;    
                    }                    
                }
		        double gamma = 1 + 1.0/planet->n; 
                if( planet->M_core >0){
                    Etot_th = - pow(gamma-1,2)/(5*M_PI/16*gamma*(3-2*gamma))*G_Newt*planet->M_core*(planet->M-planet->M_core)/planet->R_core*pow(r_RCB/planet->R_core, -(3*gamma-4)/(gamma-1));
                    //Etot_th = 8*M_PI*M_PI *G_Newt/((2-planet->n)*(2-planet->n)*(1-planet->n)) *pow(rho_RCB,2) * pow(planet->nabAd(0)*r_B, 2*planet->n)*pow(r_RCB, 2-planet->n)*pow(planet->R_core, 3- planet->n) ; 

                }
                else{
		            Etot_th = 4*M_PI*planet->Psurf[0]*(planet->scale(1),3)*(planet->n-1)/(5-planet->n) + (3.0-planet->n)/(5-planet->n) * (-G_Newt*M_Earth*M_Earth/planet->scale(1) + (planet->n+1)*planet->M*boltz_R/planet->mu_m*planet->Tsurf[0]);
                }
                fout << Etot_th  <<" " ; 
                //<< Eg << " " << Etherm << " " << Etherm2 << " " << Etot<< " " 
                fout  << m_CRB << " " << m_RCB << " " << rho_RCB << " " << r_CRB  << " " ;
                
            }
                fout << r_RCB << " " << r_solup << " " << r_liq << " " << r_crit << " " << r_soldown << " " << r_001up <<  " " << r_001down <<  " " << Eg << " " << Etherm << " " << Etot_calc << " ";
                if(planet->BL != '0'){
                    fout << planet->T_BL << " " << planet->P_BL << " ";
                }
            
                    fout << planet->M << " " << planet->M_total << " " << planet->M_dot << " " << planet->z_edge << " ";
                if(planet->R_core > 0){
                    fout << planet->M_core << " "<< planet->R_core << " "<< planet->L_core << " " ;
                }else if(planet->Jcore > 0){
                    fout << planet->M*planet->m(planet->Jcore) << " "<< planet->ym(planet->Jcore,1) << " "<< planet->ym(planet->Jcore,3) << " " << planet->ym(planet->Jcore,2) << " ";
                }
                if(planet->MassLossType == 'B'){
                    fout << planet->rho_fac << " ";
                }
                fout << m_solup << " " << m_soldown << " " << m_liq << " " << m_crit << " " << P_solup << " " << P_soldown<< " " << P_liq << " " << P_crit << " " << planet->phi(J-1) << " " << T1Gpa << " ";
                fout << planet->dt/yrs_to_s << " " << t << " " << J;
            fout.close();   
        }
        
        void set_files(struct P_P *planet ,  char num , double &f_time , std::string &_folder ,  std::string EoStype , double &Linit){
            /// Set the names of files for the particular case chosen
            char Tstring[100] , Pstring[100] , Resstring[100] , condstring[100] , corestring[100] , nstring[100] , massstring[100];
            int a;

            folder = _folder;
            std::cout << "EoS type: " << EoStype << std::endl;
            if(EoStype == "Ideal" || EoStype == "Ideal_Tab"){
                a = snprintf(nstring, 100 , "n%.1f_%.3g" , planet->n , planet->mu_m) ;
                a = snprintf(massstring, 100 , "_" );  

            }else{

                if(planet->RadioType == 'c'){
                    a = snprintf(nstring, 100 ,"MgSiO3_H%.2g_", planet->H(0));
                }else if(planet->RadioType == 'E'){
                    a = snprintf(nstring, 100 ,"MgSiO3_E_");
                }
                Minit = planet->M_total/M_Earth;
                a = snprintf(massstring, 100 , "M_E%.3g_" , Minit);  
            
            }
            if(f_time == 0){
                a = snprintf(Resstring, 100 , "Res_%.1d" , J );
            }else{
                a = snprintf(Resstring, 100 , "Res_%.1d_ftime_%.2g" , J , f_time);
            }
                
                if(planet->Tbound == 'b'){ 
                    a = snprintf(Tstring, 100 , "L%.2g_", planet->scale(3) ) ;
                }else 
                if(planet->Tbound == 'f'){
                    a = snprintf(Tstring, 100 , "Ts%.4g" , planet->Tsurf[0] );  
                    if(planet->nabtype == 'i' || planet->nabtype == 'a' ){
                        //do nothing
                    }else{
                        a += snprintf(Tstring+a, 100-a , "_L%.2g" , planet->scale(3) );
                    }
                    if(planet->MassLossType == 'B' || planet->MassLossType == 'K' ){
                        a = snprintf(Tstring+a, 100-a , "_Teq%.4g" , planet->Teq );
                    }
                }else
                if(planet->Tbound == 'i'  ){
                    a = snprintf(Tstring, 100 , "L%.2g_Teq%.4g" , planet->scale(3) , planet->Teq); 
                }else 
                if(planet->Tbound == 'm' ){
                    a = snprintf(Tstring, 100 , "L%.4g_Teq%.4gtau%.3g" , Linit , planet->Teq , planet->tau);    
                    if(planet->BL == '1'){
                        a = snprintf(Tstring+a , 100-a , "_BL1_");
                    }else 
                    if(planet->BL == '2'){
                        a = snprintf(Tstring+a , 100-a , "_BL2_");
                    }
                    else if(planet->BL == '3'){
                        a = snprintf(Tstring+a , 100-a , "_BL_");
                    }
                }else 
                if(planet->Tbound == 't' ){
                    a = snprintf(Tstring, 100 , "L%.4g_Teq%.4gtau%.3gTAB" , Linit , planet->Teq , planet->tau);   
                }else 
                if(planet->Tbound == 'l' ){
                    a = snprintf(Tstring, 100 , "L%.4g_Teq%.4gtau%.3gFIT" , Linit , planet->Teq , planet->tau);   
                }else{
                    throw std::invalid_argument("No file format found for T boundary" );
                }
                
                if(planet->Pbound == 'o'){
                    a = snprintf(Pstring, 100 , "tau%.3g",  planet->tau) ;
                }else
                if(planet->Pbound == 'c'){
                    a = snprintf(Pstring, 100 , "PRscaler%.2g", planet->opac) ;
                }else
                if(planet->Pbound == 'f'){
                    a = snprintf(Pstring, 100 , "Ps%.2g" , planet->Psurf[0] );  
                }
                a = 0;
                if(planet->M_core != 0 && planet->Jcore == 0){
                    a = snprintf(corestring, 100 , "Mc%.3g_Rc%.3g" , planet->M_core , planet->R_core );  
                    if( EoStype == "Ideal"){
                        a = snprintf(corestring+a, 100-a , "_Lc%.3g" , planet->L_core) + a;  
                        
                    }
                }
                if(planet->Jcore > 0){
                    a = snprintf(corestring+a, 100-a , "corefrac%.2g" , planet->m(planet->Jcore) );  

                }
                if(planet->MassLossType == 'c'){
                    a = snprintf(corestring+a , 100-a , "constM_dot%.2g" , planet->M_dot);
                }else if(planet->MassLossType == 'i'){
                    a = snprintf(corestring+a , 100-a , "_ISO2_M_dot_" );
                }else if(planet->MassLossType == 'P'){
                    a = snprintf(corestring+a , 100-a , "_PBC_power_" );
                }else if(planet->MassLossType == 'p'){
                    a = snprintf(corestring+a , 100-a , "_PBC_M_dot_" );
                }else if(planet->MassLossType == 'B'){
                    a = snprintf(corestring+a , 100-a , "_Booth_M_dot_" );
                }else if(planet->MassLossType == 'K'){
                    a = snprintf(corestring+a , 100-a , "_Kang_M_dot_" );
                }
                else{
                    a = snprintf(corestring+a, 100-a , "_" );
                }
            if(planet->nabtype == 'a'){
                a = snprintf(condstring, 100 , "AD" );
            }else if(planet->nabtype == 'i'){
                a = snprintf(condstring, 100 , "ISO" );
            }else
            if (planet->conds[1] == 0 && planet->conds[2] == 0){
                a = snprintf(condstring, 100 , "cond_%.3g" , planet->conds[0]);
                if(planet->ML == 'n'){
                    a = snprintf(condstring+a, 100-a , "_ML_NB" ); 
                }else if(planet->ML == 'c'){
                    a = snprintf(condstring+a, 100-a , "_ML_C" ); 
                } 
            }else{
                a = snprintf(condstring, 100 , "k0_%.3g_al%.2f_bt%.2f" , planet->conds[0], planet->conds[1] , planet->conds[2]);
            }
           
            
            bulkfile = folder+"/bulk"+num+nstring+massstring+Resstring+corestring+condstring+Tstring+Pstring+".txt";
            massfile = folder+"/massloss"+num+nstring+massstring+Resstring+corestring+condstring+Tstring+Pstring+".txt";
            
            basicfile = std::string(nstring)+massstring+Resstring+corestring+condstring+Tstring+Pstring;
            std::cout << "files set" << std::endl;

        }
        void Not_Converged(struct P_P *planet ,  EoS *eos , Converger &con  ){
            /// If the convergence algorithm has failed, abort the code and save the final values, so they can be used to figure out why
            if( eos->EoStype == "Mix"){
                Allthing = Eigen::ArrayXXd::Zero(J,I + I + 5);
                Allthing.block(0, I+3,J,1) = planet->phi;

                Allthing.block(0,I+4, J,I) = con.ErrorArray ; 
                Allthing.col(2*I+4) = ArrayXd::LinSpaced(J,0,J-1);
            }
            else{
                Allthing = Eigen::ArrayXXd::Zero(J,I + I + 3);
                
                Allthing.block(0,I+3, J,I) = con.ErrorArray ; 
            }
            Allthing.col(0) = planet->m;
            Allthing.block(0,1,J,I) = planet->y;
            Allthing.block(0, I+1,J,1) = planet->rho;
            Allthing.block(0,I+2,J,1) = planet->nabla;

            std::cout << Allthing << std::endl;
            
            int t_broken = 10;
            double time_broken = 0;
            save_bulk( planet , eos , t_broken , time_broken );
            throw std::invalid_argument("Too high an error!" );

        }
        void L_contributions(struct P_P *planet , char *timestring){
            /// Calculates various contributions to the lumminosty (for debugging)
            std::ofstream fout;
            
                Eigen::ArrayXXd Lcontribs(J,5);
                Lcontribs.col(0) = ((planet->y.col(0)*planet->scale(0)-planet->P0)/planet->P0) ;
                Lcontribs.col(1) = ((planet->y.col(2)*planet->scale(2)-planet->T0)/planet->T0);
                Lcontribs.col(2) = planet->dy.col(3);
                Lcontribs.col(3) = (planet->M /planet->scale(3) * planet->dm *  planet->H);
                Lcontribs.col(4) = (planet->M /planet->scale(3) * planet->dm *( - planet->Cp * (planet->scale(2) * planet->y.col(2)-planet->T0)/planet->dt + planet->del_rho * (planet->scale(0)*planet->y.col(0)-planet->P0)/planet->dt ));

            
        }

        Saver(int &_I, int &_J):
            I(_I) , J(_J)
        {       }
        void resize(int &_J){
            J = _J;
        }
    private:
        Eigen::ArrayXXd Fls , Allthing , others , debug_quantities , composition_array;
        std::string basicfile, folder ;
        int I , J ;
        double r_RCB =0, r_CRB = 0 , r_solup =0, r_soldown =0, r_liq = 0 , r_crit = 0 , r_001up = 0 , r_001down = 0 , m_solup = 0, m_soldown = 0, m_liq = 0 , m_crit = 0 , P_solup = 0, P_soldown = 0, P_liq = 0 , P_crit = 0 , Etot_th , Etot_calc , Eg , Etherm , Enext , Enext2 , Enext3 , Minit , T1Gpa ;
};