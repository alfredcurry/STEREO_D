#ifndef _EXTRA_MASS_TERM_H
#define _EXTRA_MASS_TERM_H

#include <iostream>
#include <Eigen/Dense>
#include "PhysicalConsts.h"
#include "StructStruct.hpp"

/// @brief Function calculating the extra heating term due to mass loss shifting mass thus changeing the potential. See MESA paper on it. 
Eigen::ArrayXd Mass_Term( P_P *planet , Eigen::ArrayXd m0 ,  Eigen::ArrayXd r0 , Eigen::ArrayXd dm0 , double M0, int j_fixed , int J){
    Eigen::ArrayXd F_m , Eg , Eg0 , E_mass_loss(planet->H.size());
    double epsilon , dm_p;
    F_m = - (planet->M*planet->m - M0*m0) / planet->dt;
    Eg = -G_Newt * planet->M * planet->m.tail(J-1) / (planet->scale(1)*planet->y.block(0,1,J-1,1));
    Eg0 = -G_Newt * M0 * m0.tail(J-1) / (r0.tail(J-1));


    for( int j=0 ; j< J-1 ; j++){
        if( j < j_fixed){
            E_mass_loss(j) = 0;
        }else{
            epsilon = (Eg(j) - Eg0(j))/planet->dt + Eg0(j)/planet->dt * ( 1 - M0*dm0(j)/(planet->M*planet->dm(j))) - (F_m(j)*Eg0(j-1)-F_m(j+1)*Eg0(j))/(planet->M*planet->dm(j));
            dm_p = M0*dm0(j) + F_m(j+1)*planet->dt + F_m(j)*planet->dt;
            E_mass_loss(j) = - epsilon /dm_p * (planet->M*planet->m(j+1)-M0*m0(j+1));
        }
    }
    return E_mass_loss;
}   

#endif