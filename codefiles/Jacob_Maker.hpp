#ifndef _JACOB_MAKER_H
#define _JACOB_MAKER_H

#include <iostream>
#include <Eigen/Dense>
#include "StructStruct.hpp"
#include "calc_nabla.hpp"
#include "EquationOfState.hpp"
#include "MeltFraction.hpp"
#include "BoundaryCondition.hpp"

class JacobMaker 
{
    //internal matrix for setting up G
    typedef Eigen::Map<Eigen::Matrix<double, 
                      Eigen::Dynamic, Eigen::Dynamic,
                      Eigen::RowMajor> > mat_t;

    public:

        double Gvec(struct P_P *planet , EoS *eos , double *Gv ); /// Function that calculates all the material values at each cell for a solution and thus finds the errors of the input solution
        void GJacob(struct P_P *planet , EoS *eos , double *l , double *d , double *u ); /// Function that calculates the jacobian (for the given solution)
    
        static Eigen::ArrayXd diffminus(const Eigen::ArrayXd x); /// Takes an array and produces an array of the same size made from the differences between array values but the 0th value takes the 0th value of the input array
        static Eigen::ArrayXd diffplus(const Eigen::ArrayXd x); /// Takes an array and produces an array of the same size made from the differences between array values but the last value takes -ve the last value of the input array
        static Eigen::ArrayXd diff(const Eigen::ArrayXd x); /// Takes an array and produces an array with one fewer coefficients made from the differences between array values

        static void prepm(struct P_P *planet ); /// Prepares the mass grid
        static void prep_material( struct P_P *planet , EoS *eos ); /// Prepares various material values, including the equation of state
        
        static Eigen::ArrayXd massbins(int J); /// Creates an array with different sized bins, in principle for an uneven mass grid
        static Eigen::ArrayXd power_law_m(int J , double X , double power); /// Creates an array where the sizes of 'cells' are determined be a power law
        static Eigen::ArrayXd power_law_Plus(int J , double X , double power , double mid); /// Creates an array where the sizes of 'cells' are determined be a power law section and a not power law section
        void resize(int _J){
            J = _J;
            nabcalc.resize(J);
        }
        B_C boundary;

        JacobMaker(int _I, int _J);

    private:
        int J , I;
        nabla_calc nabcalc;

        /// @brief Function for calculating the derivatives of the simultaneous equations that make up the jacobian
        Eigen::MatrixXd Gdiff( const Eigen::ArrayXd y  , const Eigen::ArrayXd yplus ,  const Eigen::ArrayXd scale , double M , double P0,  double T0, const double dm, const double dmm, const double rho , const double dt, const double del , const double Cp , const double nabla , double GM_R4 , double diffGM_R4 ,  double Bavg , Eigen::ArrayXXd diffnabla, double diffCpP , double diffCpT , double diffrhoP , double diffrhoT , double diffdel_rhoP, double diffdel_rhoT);
        void preprest( struct P_P *planet , EoS *eos ); /// Calculates most aspects of the set of simulataneous equations
        void prepderivs(struct P_P *planet , EoS *eos ); /// Calculates derivatives of terms for use in Gdiff
        
        //From my "pseudo atmosphere" - now not using
        double surface_factor(struct P_P *planet);
        double diffsurface_factorP(struct P_P *planet);
        
};

#endif //_JACOB_MAKER_H