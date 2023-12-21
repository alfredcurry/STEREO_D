#ifndef _SURFACE_PRESS_H
#define _SURFACE_PRESS_H

#include <iostream>
#include <Eigen/Dense>
#include "PhysicalConsts.h"

template <class I>
class P_surf{
    public:
        static I press(const I M , const I R , const I T , const I kap , const double al_opac , const double bt_opac , const double tau )
        {   
            return  pow( tau * G_Newt * M* (al_opac + 1) /(R*R*kap * pow(T, bt_opac)) , 1.0/(al_opac+1));  
        }
        static I diff_R(const I Ps, const I R , const double al_opac)
        {
            return -  2 / (al_opac + 1) * Ps/R;
        }
        static I diff_Teff(const I Ps, const I T , const double al_opac , const double bt_opac )
        {
            return - bt_opac / (al_opac + 1) * Ps/T;
        }
};
#endif
