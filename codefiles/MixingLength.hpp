#ifndef _MIXING_l_H
#define _MIXING_l_H

template <typename Real> /// Class for calculatin the mixing length via different methods
class MixingLength
{
    private:
            /* data */
    public:
        static Real Nearest_Bound(double R , Real r){
            return std::min(r , R - r); 
        }
        static Real Scale_Height(Real P , Real g , Real rho){
            return P/(g*rho);
        }
        static Real Constant(double R){
            return R/4;
        }

};             
        

#endif