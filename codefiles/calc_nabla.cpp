#include <Eigen/Dense>
#include "adiabatnab.hpp"
#include "Jacob_Maker.hpp"
#include "StructStruct.hpp"
#include "Boost_solver.hpp"
#include "ConductiveNab.hpp"
#include "Conductivity.hpp"
#include "calc_nabla.hpp"
#include "EquationOfState.hpp"
#include "MixingLength.hpp"
#include "Boundary_Layer.hpp"

using namespace Eigen; 

void nabla_calc::nabla( struct P_P *planet){
    if(planet->nabtype != 'i' && planet->nabtype != 'a' && planet->nabtype != 'c' && planet->nabtype != 'B' && planet->nabtype != 'g' ){
        std::cout << "nabtype: " << planet->nabtype << std::endl;
        throw std::invalid_argument("Invalid nabtype. Valid types: 'i', 'a', 'c', 'B' or 'g'. ");
    }
    if (planet->nabtype == 'i'){
        planet->nabla = ArrayXd::Zero(planet->nabla.size());    
    }else{
        planet->nabAd = Nab_Ad<ArrayXd>::nab(planet->ym.col(0) , planet->ym.col(2) , planet->del_rho , planet->Cp );
        if (planet->nabtype == 'a'){
            planet->nabla = planet->nabAd;
        }else{ 
            Re_crit = 9.0/8;
		
            double A , B , C;     
            if(planet->Jcore > 0){
                planet->nabla.head(planet->Jcore) = planet->nabAd.head(planet->Jcore);
                planet->supernabla.head(planet->Jcore) = 0;
                planet->nabNonAd.head(planet->Jcore) = 0;
                planet->L_cond.head(planet->Jcore) = 0;
                planet->L_conv.head(planet->Jcore) = 0;
                planet->L_grav.head(planet->Jcore) = 0;
                planet->L_mix.head(planet->Jcore) = 0;
                planet->u_conv.head(planet->Jcore) = 0;
                planet->l_mix.head(planet->Jcore) = 0;
            }
            for(int j=planet->Jcore ; j<J-1 ; j++){
                //if (planet->ym(j,3) == 0 || planet->nabtype == 'n'){
                  //  planet->nabla(j) = planet->nabAd(j);
                    //planet->L_conv(j) = planet->ym(j,3);
                    //planet->L_cond(j) = 0;
                //}

                if(planet->nabtype == 'c' || planet->phi(j) >= 1 || planet->phi(j) <= 0){
                    planet->L_cond(j) = planet->ym(j,3);
                    planet->L_grav(j) = 0;
                }else{
                    C = C_grav(planet , j);
                    planet->L_grav(j) = Fluxes<double>::LGrav(C, planet->phi(j), planet->DeltaS(j));
                    planet->L_cond(j) = planet->ym(j,3) - planet->L_grav(j);
                }

                planet->nabla(j) = CondNab<double>::nab(planet->L_cond(j) , planet->K(j));
                planet->supernabla(j) = planet->nabla(j) - planet->nabAd(j);

                if (planet->supernabla(j)>0){
                    
                    if (planet->K(j) == 0){
                        planet->supernabla(j) = 100;
                    }
                    if(planet->ML == 'c'){
                        planet->l_mix(j) = MixingLength<double>::Constant(planet->scale(1));
                    }else if (planet->ML == 'n'){
                        planet->l_mix(j) = std::max(MixingLength<double>::Nearest_Bound(planet->R0+planet->z_edge-planet->R_core0, planet->ym(j,1)-planet->R_core0) , 2e3);
                    }
                    
                    if (planet->l_mix(j) <= 0){
                        
                        std::cout << j <<  " MIXING LENGTH PROBLEM " << planet->l_mix(j) << std::endl;
                    }
    
                    AB( planet , A , B  , j);
                    planet->nabNonAd(j) = nabnonadcalc(planet , A , B , C , Re_crit ,  j);

                    if(planet->nabtype == 'c' || planet->phi(j) >= 1 || planet->phi(j) <= 0){
                        Fluxes<double> Flux(Re_crit , A , B , planet->l_mix(j)/planet->nu(j), planet->K(j) , planet->nabAd(j));
                        
                        planet->L_conv(j) = Flux.LConv(planet->nabNonAd(j));
                        planet->L_cond(j) = Flux.LCond(planet->nabNonAd(j));
                        planet->L_grav(j) = 0;
                        planet->L_mix(j) = 0;
                        planet->u_conv(j) = Flux.u ;
                    
                    }else{    
                        Fluxes<double> Flux(Re_crit , A , B , planet->l_mix(j)/planet->nu(j), planet->K(j) , planet->nabAd(j), planet->Cp(j) , planet->ym(j,2)*planet->diffphiT(j) , planet->ym(j,0)*planet->diffphiP(j) );

                        planet->L_conv(j) = Flux.LConv(planet->nabNonAd(j));
                        planet->L_cond(j) = Flux.LCond(planet->nabNonAd(j));
                        planet->L_mix(j) = 0;
                        planet->u_conv(j) = Flux.u ;
                        
                        if( planet->nabtype == 'B'){
                            planet->L_mix(j) = Flux.LMix(planet->nabNonAd(j) , planet->DeltaS(j));
                        }    
                    }
                    
                    planet->nabla(j) = planet->nabAd(j) + planet->nabNonAd(j);
                    
                    
                }else{
                    
                    planet->L_cond(j) = planet->ym(j,3);
                    planet->L_conv(j) = 0;
                    planet->L_grav(j) = 0;
                    planet->L_mix(j) = 0;
                    planet->nabNonAd(j) = 0;
                    planet->u_conv(j) = 0;
                    planet->l_mix(j) = 0;
                }
            }
        }}

    }
void nabla_calc::diffnabla( struct P_P *planet0 , EoS *eos){
    if (planet0->nabtype == 'i'){
        planet0->diffnabla =0;
    }else{
    
        if(planet0->nabtype == 'a'){
            planet0->diffnabla.col(0) = Nab_Ad<ArrayXd>::DiffnabP(planet0->ym.col(0) , planet0->ym.col(2) , planet0->del_rho , planet0->Cp ,  planet0->diffdel_rhoP, planet0->diffCpP) ;
            planet0->diffnabla.col(1) = 0;
            planet0->diffnabla.col(2) = Nab_Ad<ArrayXd>::DiffnabT(planet0->ym.col(0) , planet0->ym.col(2) , planet0->del_rho , planet0->Cp ,  planet0->diffdel_rhoT, planet0->diffCpT);
            planet0->diffnabla.col(3) = 0;

        }else{
            CondNab<double> cond_nab;
            if(planet0->Jcore > 0){
                planet0->diffnabla.block(0,0,planet0->Jcore,1) = Nab_Ad<ArrayXd>::DiffnabP(planet0->ym.block(0,0,planet0->Jcore,1) , planet0->ym.block(0,2,planet0->Jcore,1) , planet0->del_rho.head(planet0->Jcore) , planet0->Cp.head(planet0->Jcore) ,  planet0->diffdel_rhoP.head(planet0->Jcore), planet0->diffCpP.head(planet0->Jcore));
                planet0->diffnabla.block(0,1,planet0->Jcore,1) = 0;
                planet0->diffnabla.block(0,2,planet0->Jcore,1) = Nab_Ad<ArrayXd>::DiffnabT(planet0->ym.block(0,0,planet0->Jcore,1) , planet0->ym.block(0,2,planet0->Jcore,1) , planet0->del_rho.head(planet0->Jcore) , planet0->Cp.head(planet0->Jcore) ,  planet0->diffdel_rhoT.head(planet0->Jcore), planet0->diffCpT.head(planet0->Jcore));
                planet0->diffnabla.block(0,3,planet0->Jcore,1) = 0;
            }
            for(int j=planet0->Jcore ; j<J ; j++){
                
                if (planet0->supernabla(j)<=0 && (planet0->nabtype == 'c' || planet0->phi(j) >= 1 || planet0->phi(j) <= 0) ){
                    
                    planet0->diffnabla(j,0) =  cond_nab.diffnabP(planet0->ym(j,3) , planet0->K(j) , planet0->ym(j,0) , planet0->rho(j) , planet0->diffrhoP(j), planet0->k_cond(j), Conduct<double>::DiffkP(planet0->rho(j), planet0->diffrhoP(j) ,  planet0->ym(j,2) , planet0->conds[0], planet0->conds[1], planet0->conds[2]));
                    planet0->diffnabla(j,1) =  0;    
                    planet0->diffnabla(j,2) =  cond_nab.diffnabT(planet0->ym(j,3) , planet0->ym(j,2) , planet0->K(j) , planet0->rho(j) , planet0->diffrhoT(j), planet0->k_cond(j), Conduct<double>::DiffkT(planet0->rho(j), planet0->diffrhoT(j) , planet0->ym(j,2) , planet0->conds[0], planet0->conds[1], planet0->conds[2]));
                    planet0->diffnabla(j,3) =  cond_nab.diffnabL(planet0->K(j));    
                    
                }else{

                    ArrayXd dy = 1e-10 * planet0->ym.row(j);
                    P_P planet(I,1);
                    planet.scale = planet0->scale;
                    planet.M = planet0->M;
                    planet.mu_m = planet0->mu_m;
                    planet.n = planet0->n;
                    planet.nabtype = planet0->nabtype;
                    planet.a_crystal = planet0->a_crystal;
                    planet.P0(0) = planet0->P0(j);
                    planet.T0(0) = planet0->T0(j);
                    planet.Jcore = 0;
                    
                    planet.mm(0) = planet0->mm(j);
                    
                    for (int i=0; i<3 ; i++){ 
                        planet.conds[i] = planet0->conds[i];
                    }
                    if(j==J-1){
                        planet.y(0,1) = planet0->ym(j-1,1);
                    }
                    double  A , B , C;

                    for(int i=0; i < I ; i++){ 
                        
                        planet.ym = planet0->ym.row(j);
                        
                        planet.ym(i) = planet.ym(i) + dy(i);

                        JacobMaker::prep_material(&planet , eos);
                        C = C_grav(&planet , 0);

                        if(planet0->nabtype == 'c' || planet.phi(0) >= 1 || planet.phi(0) <= 0){
                            planet.L_cond(0) = planet.ym(3);
                        }else{
                            planet.L_grav(0) = Fluxes<double>::LGrav(C, planet.phi(0), planet.DeltaS(0));
                            planet.L_cond(0) = planet.ym(3) - planet.L_grav(0);
                        }
                        planet.nabla = CondNab<double>::nab(planet.L_cond(0) , planet.K(0));
                        planet.nabAd = Nab_Ad<double>::nab(planet.ym(0) , planet.ym(2) , planet.del_rho(0) , planet.Cp(0) );
                        
                        planet.supernabla = planet.nabla - planet.nabAd;
                        if( planet.supernabla(0) <= 0 ){
                            std::cout << j << " SWITCHED BELOW ADIABAT " << i << std::endl;
                            planet0->diffnabla(j,i) = (planet.nabla(0) - planet0->nabla(j))/dy(i);
                        }
                        else{
                            if (planet.K(0) == 0){
                                planet.supernabla(0) = 100;
                            }    

                            if(planet0->ML == 'c'){
                                planet.l_mix = MixingLength<double>::Constant(planet0->scale(1));
                            }else if (planet0->ML == 'n'){
                                planet.l_mix = std::max(MixingLength<double>::Nearest_Bound(planet0->R0+planet0->z_edge-planet0->R_core0, planet.ym(1)-planet0->R_core0) , 2e3);
                            }
                            if (planet.l_mix(0) <= 0){
                        
                                std::cout <<  " MIXING LENGTH PROBLEM " << planet.l_mix(0) << std::endl;
                            }
                            //planet.l_mix = MixingLength<double>::Constant(planet.scale(1));
                            AB( &planet , A , B , 0);
                            planet.nabNonAd = nabnonadcalc(&planet , A , B , C , Re_crit , 0);                            
                            
                            planet0->diffnabla(j,i) = (planet.nabNonAd(0) - planet0->nabNonAd(j))/dy(i);
                            
                            if(i==0){                                
                                planet0->diffnabla(j,i) = planet0->diffnabla(j,i) +  Nab_Ad<double>::DiffnabP(planet0->ym(j,0) , planet0->ym(j,2) , planet0->del_rho(j) , planet0->Cp(j) , planet0->diffdel_rhoP(j) , planet0->diffCpP(j));
                            }
                            else if(i==2){
                                planet0->diffnabla(j,i) = planet0->diffnabla(j,i) + Nab_Ad<double>::DiffnabT(planet0->ym(j,0) , planet0->ym(j,2) , planet0->del_rho(j) , planet0->Cp(j) , planet0->diffdel_rhoT(j) , planet0->diffCpT(j));
                            }

                        }

                    }
                }
            }   
        }
    }

            planet0->diffnabla.col(0) = planet0->scale(0) * planet0->diffnabla.col(0) ;
            planet0->diffnabla.col(1) = planet0->scale(1) * planet0->diffnabla.col(1) ;
            planet0->diffnabla.col(2) = planet0->scale(2) * planet0->diffnabla.col(2) ;
            planet0->diffnabla.col(3) = planet0->scale(3) * planet0->diffnabla.col(3) ;

            //std::cout <<"\n P T\n" <<  planet0->ym<< std::endl;
            //std::cout <<"diffnabla\n" <<  planet0->diffnabla<< std::endl;
}    


nabla_calc::nabla_calc(int _I , int _J)
        :   
            I(_I) , J(_J)
        {};
void nabla_calc::AB(struct P_P *planet , double &A , double &B , int j){
    A = planet->del_rho(j)*pow(planet->g(j),2)*planet->rho(j)*planet->rho(j)*pow(planet->l_mix(j),3)/(18*planet->nu(j)*planet->ym(j,0));
    B = 4*M_PI*G_Newt*planet->M*planet->mm(j)*pow(planet->rho(j),2)*planet->Cp(j)*planet->ym(j,2)*planet->l_mix(j)/planet->ym(j,0); 
    
}
double nabla_calc::C_grav(struct P_P *planet , int j  ){
    return -4*M_PI*G_Newt*planet->M*planet->mm(j)*planet->a_crystal*planet->a_crystal * planet->rho(j)*(planet->rho_m(j)-planet->rho_s(j))*planet->ym(j,2)/planet->visc_m;
}
double nabla_calc::nabnonadcalc(struct P_P *planet , const double &A , const double &B , const double &C , double Re_crit , int j ){
    /// Uses Algorithm TOMS 748 (essentially Brent's method) to find the super-adiabatic part of the logTlogP gradient
            
            if (planet->supernabla(j) < 0){
                std::cout <<"oh dear: "<< planet->supernabla << std::endl;
            }
            if(planet->nabtype == 'c' || planet->phi(j) >= 1 || planet->phi(j) <= 0){
                //std::cout << j << " a problem to TOM " <<  A << " " << planet->nu(j) << " " << planet->rho(j) << std::endl;
                return boost_Tom(planet->ym(j,3) , Re_crit , A , B , planet->l_mix(j)/planet->nu(j) , planet->K(j) ,  planet->nabAd(j) , planet->supernabla(j));   
            }
            else{
                //std::cout <<j << " " <<  A << " " << B << " " << planet->nabAd(j) <<" " << planet->K(j) << " " << C << " " << planet->DeltaS(j) << " " << planet->diffphiT(j) << std::endl;
                return boost_Tom(planet->ym(j,3) , Re_crit , A , B , planet->l_mix(j)/planet->nu(j) , planet->K(j) ,  planet->nabAd(j) , planet->DeltaS(j) , planet->phi(j), planet->Cp(j) , planet->ym(j,2)*planet->diffphiT(j) , planet->ym(j,0)*planet->diffphiP(j) , C , 1.0 , planet->nabtype);   
            }
        }

