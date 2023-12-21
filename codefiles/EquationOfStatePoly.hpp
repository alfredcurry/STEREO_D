#ifndef _EQUATION_OF_STATE_H
#define _EQUATION_OF_STATE_H

#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

void EoSrhopoly(const Ref<ArrayXd> P, ArrayXd &rho , double K, double n)
{   
    rho =  pow(P / K , n/(n+1)) ;  
}

void DiffrhoPpoly(const Ref<ArrayXd> P, const Ref<ArrayXd>rho ,  ArrayXd &diffrho , double n)
{ 
    diffrho = n/(n+1) * rho/P ;
}
void DiffrhoTpoly(const Ref<ArrayXd> T, const Ref<ArrayXd>rho ,  ArrayXd &diffrho , double n)
{ 
    diffrho = n* rho/T ;
}
void Diff2rhoTTpoly(const Ref<ArrayXd> T, const Ref<ArrayXd>rho ,  ArrayXd &diffrho , double n)
{ 
    diffrho = n*(n-1) * rho/(T*T) ;
}
void Diff2rhoPTpoly(const Ref<ArrayXd> P, const Ref<ArrayXd> T, const Ref<ArrayXd>rho ,  ArrayXd &diffrho , double n)
{ 
    diffrho = -n/(n+1) * rho/(P*T) ;
}

#endif 