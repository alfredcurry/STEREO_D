#ifndef _VOLCANIC_H
#define _VOLCANIC_H

#include <Eigen/Dense>

template <typename Real>
Real Volcanic_Luminosity( Eigen::Array<Real , Eigen::Dynamic , 1 > &phi , Eigen::Array<Real , Eigen::Dynamic , 1 > &u , Eigen::Array<Real , Eigen::Dynamic , 1 > &dm ){
    return (phi*u*dm).sum();
}

#endif