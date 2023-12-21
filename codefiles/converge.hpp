#ifndef _CONVERGE_H
#define _CONVERGE_H

#include <iostream>
#include <Eigen/Dense>
#include "StructStruct.hpp"
#include "Jacob_Maker.hpp"
#include "EquationOfState.hpp"

class Converger
{
    /// Class for the convergence algorithm for the Henyey scheme
    public:
    
        double *Gv , maxinacc , meaninacc;
        int n , nmax;
        Eigen::ArrayXXd ErrorArray;

        /// Function that runs the convergence procedure
        void converge(struct P_P *planet , EoS *eos , double acc);
        Converger(int _I , int _J);
        void resize(int _J){
            /// Resizes various arrays when the grid-size changes
            J = _J;
            Jacob.resize(J);
            planetplus.resize(J);
            delete x;
            delete l;
            delete d;
            delete u;
            delete Gv;
            delete Gvplus;
            x = new double[I*J] ;
            l = new double[I*I*J]; 
            d = new double[I*I*J] ;
            u = new double[I*I*J];
            Gv = new double[I*J];
            Gvplus = new double[I*J];
        }
        JacobMaker Jacob;

    private:
        int I , J, again;
        double *x, *l, *d , *u , *Gvplus;
        P_P planetplus; 
        
        /// Adds the computed correction to the solution
        void add_correction(struct P_P *planet , Eigen::ArrayXXd Dy );
        /// Function for checking the analytically calculated derivatives in the Jacobian vs. brute force calculated versions
        void checker(struct P_P *planet , EoS *eos , Eigen::MatrixXd G , Eigen::MatrixXd diffG , double acc , int row , int Gtype , int col );
        /// Calculates the magnitude of the gradient from the Jacobian in the direction of the correction
        double gradient(double *x, const Eigen::MatrixXd &diffGd , const Eigen::MatrixXd &diffGu , const Eigen::MatrixXd &diffGl , const Eigen::MatrixXd &G);
        /// Searches for the lowest error solution along the direction of the correction
        bool lsrch(struct P_P *planetplus , EoS *eos , double &f , double *Gv , const double slope , const Eigen::ArrayXXd Dy , const struct P_P *planet );

};

#endif //_CONVERGE_H