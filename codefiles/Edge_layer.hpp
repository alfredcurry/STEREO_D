#ifndef _EDGE_H
#define _EDGE_H

#include <vector>
#include "PhysicalConsts.h"

template <class Real>
Real Edge_depth( const Real &z_min ,  const Real &frac , const Real &T0 , const Real &Tss , const Real &k_cond){
    return std::max(2*k_cond*(4.0/5*Tss - T0 )/(frac*Sig_Stefan*pow(Tss,4)), z_min);
}
template <class Real>
Real Edge_pressure( const Real &z_depth ,  const Real &rho, const Real &g , const Real &P0_min){
    return std::min(z_depth*rho*g , P0_min);
}

#endif