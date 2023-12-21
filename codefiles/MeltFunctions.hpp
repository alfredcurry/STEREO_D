#include <Eigen/Dense>
#include "MeltFraction.hpp"
#include <boost/math/interpolators/cubic_hermite.hpp>
#include <fstream>
#include <vector>
#include "Vector_reader.h"
#include <cmath>
#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>
#include <boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_fehlberg78.hpp>

/// @brief  Implementation of MeltFraction class. Often phi is calculated simply as the function we write as NormalisedT. However, this would not have a smooth derivative. Therefore we use a smooth approximation to a top hat function, multiplied by NormalisedT, rather than NormalisedT itself.

template <typename Real>
Real NormalisedT(Real T , Real Tsol , Real Tliq){
    return (T - Tsol) / (Tliq - Tsol);
}
template <typename Real>
Real NormalisedTdPdT(Real T , Real Tsol , Real Tliq , Real diffTsoldP , Real diffTliqdP){
    return -1.0/pow(Tliq-Tsol,2) * (diffTliqdP - diffTsoldP) ;
}
template <typename Real>
Real NormalisedTdP(Real T , Real Tsol , Real Tliq , Real diffTsoldP, Real diffTliqdP){
    return -diffTsoldP /(Tliq-Tsol) + (T - Tsol)*NormalisedTdPdT(T,Tsol,Tliq,diffTsoldP,diffTliqdP);
}
template <typename Real>
Real NormalisedTdT(Real T , Real Tsol , Real Tliq){
    return 1.0 / (Tliq - Tsol);
}

/// Smooth approximation to top hat function using a gaussian- no longer used
template <typename Real>
Real SmoothTopHat(Real x  , double n , Real A){
    return exp( - n /(1.0-pow(2.0*x-1.0 , 2)) ) / A;
}
template <typename Real>
Real SmoothTopHatdx(Real x , double n , Real A){
    return -4.0*(2.0*x-1.0)*n/pow( 1.0-pow(2.0*x-1.0,2) , 2) * exp( - n /(1.0-pow(2.0*x-1.0,2)) ) / A;
}

typedef double state_type;

class TopHat_integrator {
    private:
        double n , A;
    public:
        void operator() ( const state_type &phi , state_type &dphidx , const double x ){
            dphidx = SmoothTopHat(x , n , A);
        }
        TopHat_integrator(double _n , double _A) : 
            n(_n) , A(_A) {}
};

void integrate_top_hat( state_type &phi , const double x0 , const double x1, const double n , const double A){
    // Error stepper, used to create the controlled stepper
    typedef boost::numeric::odeint::runge_kutta_fehlberg78< double > rkck54;

    // Controlled stepper: 
    // it's built on an error stepper and allows us to have the output at each 
    // internally defined (refined) timestep, via integrate_adaptive call 
    typedef boost::numeric::odeint::controlled_runge_kutta< rkck54 > ctrl_rkck54;
    
    TopHat_integrator hat(n , A);

    boost::numeric::odeint::integrate_adaptive( ctrl_rkck54() , hat , phi, x0 , x1 ,  1e-3 );
    
}


template <class Real> /// Smooth approximation to a top hat function using a polynomoal
class PolyTopHat {
    private:
        Real A , b , c , sigma;
    public:
        Real operator()(Real &x){
            if(x < sigma){
                return b*x*x + c*x*x*x;
            }else
            if(x >1-sigma){
                return b*pow(1-x,2) + c*pow(1-x,3);
            }else{
                return A;
            }
        }   
        Real integral(Real &x){
            if(x < sigma){
                return b*x*x*x/3 + c*pow(x,4)/4;
            }else
            if(x >1-sigma){
                return - b*pow(1-x,3)/3 - c*pow(1-x,4)/4 + 1;
            }else{
                return A*x - A*sigma/2;
            }
        }
        Real prime(Real &x){
            if(x < sigma){
                return 2*b*x + 3*c*x*x;
            }else
            if(x >1-sigma){
                return - 2*b*(1-x) - 3*c*pow(1-x,2);
            }else{
                return 0;
            }
        } 
        PolyTopHat(Real _sigma):
        sigma(_sigma){
            A = 1.0/(1.0-sigma);
            b = 3.0*A/(sigma*sigma);
            c = -2.0*A/(sigma*sigma*sigma);
        }
};

template <typename Real>
struct SimonGlatzel{//Simon Glatzel fit
    public:
        static Real eval(Real const &p , Real &T0 , Real &a , Real &c){
           return T0*pow((p) / a + 1.0 , 1.0/c);
        }
        static Real prime(Real const &p , Real &T0 , Real &a , Real &c){
            return T0/(a*c) * pow((p) / a + 1.0 , 1.0/c - 1);
        }
};



template <typename Real>
Eigen::Array<Real , Eigen::Dynamic , 1 > MeltFraction<Real>::Phi(Eigen::Array<Real , Eigen::Dynamic , 1 > P, Eigen::Array<Real , Eigen::Dynamic , 1 > T){

    Eigen::Array<Real , Eigen::Dynamic , 1 > phi ;
    phi = Eigen::ArrayXd::Zero(P.size()) ;
    Tsol = Eigen::ArrayXd::Zero(P.size()) ;
    Tliq = Tsol;
    NormT = 0*Tsol;
    TsoldP = 0*P ;
    TliqdP = 0*P ;
    diff2phiPT.resize(T.size());  
    diff2phiTT.resize(T.size());  
    NormTdT.resize(T.size());

    PolyTopHat<double> TopHat(smooth_parameters[0]);
    for(int j=0 ; j < P.size() ; j++){
        Tsol(j) = Solidus(P(j));
        Tliq(j) = Liquidus(P(j));
        if(Tsol(j) > Tliq(j)){/// This can happen because the fit doesn't apply to infinte P. Generally we are not close to the solidus at this P, but a value does need to be calculated
            Tsol(j) = (Tliq(j)+Tsol(j))/2;
            Tliq(j) = Tsol(j);
            if( T(j) > Tliq(j)){
                phi = 1;
                NormT = 2;
            }else{
                phi = 0;
                NormT = -1;
            }
        }else{
            TsoldP(j) = DiffTsoldP(P(j));
            TliqdP(j) = DiffTliqdP(P(j));
            NormT(j) = NormalisedT(T(j) , Tsol(j) , Tliq(j));
            if(NormT(j) <= 0){
                phi(j) = 0;
            }else
            if(NormT(j) >=1){
                phi(j) = 1;
            }
            else{
                phi(j) = TopHat.integral(NormT(j)); /// 
                if(phi(j) > 1.0){ //This may arise due to numerical accuracy (although no more because we now use a fully analytical version)
                    std::cout << NormT(j) << " integration problem " <<std::endl;
                    phi(j) = 1;
                }
            }
        }

    }
    return phi;
}
template <typename Real>
Eigen::Array<Real , Eigen::Dynamic , 1 > MeltFraction<Real>::DiffPhiT(Eigen::Array<Real , Eigen::Dynamic , 1 > T){
    Eigen::Array<Real , Eigen::Dynamic , 1 > diffPhiT = 0*T;
    //std::cout << NormTdT.size() << " " << T.size() << std::endl;

    PolyTopHat<double> TopHat(smooth_parameters[0]);
    for(int j=0 ; j < T.size() ; j++){ 
        if(NormT(j) > 0 && NormT(j)<1 ){
            NormTdT(j) = NormalisedTdT(T(j) , Tsol(j) , Tliq(j)) ;
            diffPhiT(j) = NormTdT(j) * TopHat(NormT(j));
//            diffPhiT(j) = NormTdT(j) * SmoothTopHat(NormT(j) , smooth_parameters[0] , smooth_parameters[1] );
        }   
    }
    return diffPhiT;
}
template <typename Real>
Eigen::Array<Real , Eigen::Dynamic , 1 > MeltFraction<Real>::DiffPhiP( Eigen::Array<Real , Eigen::Dynamic , 1 > P , Eigen::Array<Real , Eigen::Dynamic , 1 > T ){

    Eigen::Array<Real , Eigen::Dynamic , 1 > diffPhiP = 0*P ;
    NormTdP = 0*P;
    PolyTopHat<double> TopHat(smooth_parameters[0]);
    
    for(int j=0 ; j < P.size() ; j++){
        if(NormT(j) > 0 && NormT(j)<1 ){
            
            NormTdP(j) = NormalisedTdP(T(j) , Tsol(j) , Tliq(j) , TsoldP(j) , TliqdP(j)) ;            
            diffPhiP(j) = NormTdP(j) * TopHat(NormT(j));
//            diffPhiP(j) = NormTdP(j) * SmoothTopHat(NormT(j) , smooth_parameters[0] , smooth_parameters[1] );  
            
        }
    }
   
    return diffPhiP /1e9;
}

template <typename Real>
void MeltFraction<Real>::SecondDerivs(Eigen::Array<Real , Eigen::Dynamic , 1 > P , Eigen::Array<Real , Eigen::Dynamic , 1 > T){
    
    PolyTopHat<double> TopHat(smooth_parameters[0]);

    for(int j=0 ; j < T.size() ; j++){
        if(NormT(j) > 0 && NormT(j)<1 ){
            
            diff2phiTT(j) = TopHat.prime(NormT(j)) * NormTdT(j)*NormTdT(j) ;
            diff2phiPT(j) = TopHat.prime(NormT(j)) * NormTdT(j)*NormTdP(j) + TopHat(NormT(j))*NormalisedTdPdT(T(j),Tsol(j),Tliq(j),TsoldP(j),TliqdP(j)) ;
  
  //          diff2phiTT(j) = SmoothTopHatdx(NormT(j) , smooth_parameters[0] , smooth_parameters[1] ) * NormTdT(j)*NormTdT(j) ;
    //        diff2phiPT(j) = SmoothTopHatdx(NormT(j) , smooth_parameters[0] , smooth_parameters[1] ) * NormTdT(j)*NormTdP(j) + SmoothTopHat(NormT(j) , smooth_parameters[0] , smooth_parameters[1] )*NormalisedTdPdT(T(j),Tsol(j),Tliq(j),TsoldP(j),TliqdP(j)) ;
        
        }
    }
    diff2phiPT = diff2phiPT/1e9;
}
template <typename Real>
Real MeltFraction<Real>::Solidus( Real P){
    if(P < 100e9){
        return Solidus_spline(P/1e9);
    }else{
        return SimonGlatzel<Real>::eval(P/1e9,analytical_parameters[0] , analytical_parameters[1] , analytical_parameters[2]);
    }
}
template <typename Real>
Real MeltFraction<Real>::Liquidus( Real P){
    if(P < 100e9){
        return Liquidus_spline(P/1e9);
    }else{
        return SimonGlatzel<Real>::eval(P/1e9,analytical_parameters[3] , analytical_parameters[4] , analytical_parameters[5]);
    }
}
template <typename Real>
Real MeltFraction<Real>::DiffTsoldP(Real P){
    if(P < 100e9){
        return Solidus_spline.prime(P/1e9);
    }else{
        return SimonGlatzel<Real>::prime(P/1e9,analytical_parameters[0] , analytical_parameters[1] , analytical_parameters[2]);
    }
}
template <typename Real>
Real MeltFraction<Real>::DiffTliqdP(Real P){
    if(P < 100e9){
        return Liquidus_spline.prime(P/1e9);
    }else{
        return SimonGlatzel<Real>::prime(P/1e9,analytical_parameters[3] , analytical_parameters[4] , analytical_parameters[5]);
    }
}

template <typename Real>
Real MeltFraction<Real>::DeltaS(Real DeltaV , Real phi , Real TsoldP , Real TliqdP){
    /// @brief see Appendix of Curry et. al. (2023)
    return (phi/TliqdP + (1-phi)/TsoldP ) * DeltaV; // If constant use  39.9+-3.3 J/K/mol (Lesher & Spera, 2015) - const. not bad approx. see Stixrude (2009) Fig. 2
}
template <typename Real>
Real MeltFraction<Real>::DiffDeltaSdP(Real DeltaS , Real DeltaV , Real DeltaVdP , Real diffphidP , Real TsoldP , Real TliqdP){
    return DeltaS * DeltaVdP / DeltaV + DeltaV*diffphidP*(1.0/TliqdP - 1.0/TsoldP);
}
template <typename Real>
Real MeltFraction<Real>::DiffDeltaSdT(Real DeltaS , Real DeltaV , Real DeltaVdT , Real diffphidT , Real TsoldP , Real TliqdP){
    return DeltaS * DeltaVdT / DeltaV + DeltaV*diffphidT*(1.0/TliqdP - 1.0/TsoldP);
}

template <typename Real> /// Constructor reads in the data for the solidus/liquidus
MeltFraction<Real>::MeltFraction():
    smooth_parameters(2) , analytical_parameters(6)
{
    std::string file , filebase = "eos/SolLiq/";
    char filestring[100];
    int ROWS = 2000;
    double P0 = 0, dP = 0.1;
   
    std::vector<Real> Liq_data(ROWS) , Sol_data(ROWS) , Liq_dP_data(ROWS) , Sol_dP_data(ROWS) ;
    snprintf(filestring , 100 , "And11_Litasov_%.d_%.1f.txt" , ROWS , dP);
    
    file = filebase+"Solidus/"+filestring;
    std::cout << "Reading Liquidus..." << std::endl;
    readVector(Sol_data , ROWS , file);
    
    file = filebase+"Liquidus/"+filestring;
    std::cout << "Reading Solidus..." << std::endl;
    readVector(Liq_data , ROWS , file);

    snprintf(filestring , 100 , "And11_LitasovdP_%.d_%.1f.txt" , ROWS , dP);
    
    file = filebase+"Solidus/"+filestring;
    std::cout << "Reading Liquidus derivative..." << std::endl;
    readVector(Sol_dP_data , ROWS , file);
    
    file = filebase+"Liquidus/"+filestring;
    std::cout << "Reading Solidus..." << std::endl;
    readVector(Liq_dP_data , ROWS , file);

    Liquidus_spline = boost::math::interpolators::cardinal_cubic_hermite<std::vector<Real>>(std::move(Liq_data), std::move(Liq_dP_data), P0, dP);
    Solidus_spline = boost::math::interpolators::cardinal_cubic_hermite<std::vector<Real>>(std::move(Sol_data), std::move(Sol_dP_data), P0, dP);
    
    double Sol_T0=2045.0  , Sol_a=92.0 , Sol_c=1.3 , Liq_T0=1940.0 , Liq_a=29.0 , Liq_c=1.9;
    analytical_parameters = {Sol_T0  , Sol_a , Sol_c , Liq_T0 , Liq_a , Liq_c};

    smooth_parameters[0] = 0.1;
 //   smooth_parameters[1] = 0;

//    integrate_top_hat( smooth_parameters[1] , 0.0 , 1.0 , smooth_parameters[0] , 1.0); 
  //  std::cout << "A = " << smooth_parameters[1] << std::endl;
    std::cout << "Solidus / Liquidus interpolators complete" << std::endl;
}
