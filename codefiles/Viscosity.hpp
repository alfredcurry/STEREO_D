#ifndef _VISCOSITY_H
#define _VISCOSITY_H

#include <vector>

template <typename Real> /// Class for calculating the viscosity its derivatives w.r.t T, P and phi
class Viscosity
{
    private: 
        std::vector<Real> variables;

    public:
        Real Visc(Real &phi , Real &P , Real &T);
        Real DiffVisc(Real &phi);
        Real DiffVisc_dT(Real &phi, const Real &visc , Real &P , Real &T , Real &dphiT);
        Real DiffVisc_dP(Real &phi, const Real &visc , Real &P , Real &T , Real &dphiP , Real &dT_phi0dP);

        Viscosity(std::vector<Real> _variables);
};

#endif