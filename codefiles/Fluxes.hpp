template <class T> /// Heat flux for a super-adiabatic gradient
struct Fluxes
{
    Fluxes( T const & _Re_crit , T const & _A , T const & _B , T const & _l_nu, T const & _K_cond , T const & _nab_ad , T const & _Cp , T const & _TdiffphiT , T const & _PdiffphiP ) : 
       Re_crit(_Re_crit) , A_vis(_A) ,  B_conv(_B) , l_nu(_l_nu) , K_cond(_K_cond) , nab_ad(_nab_ad) , Cp(_Cp) , TdiffphiT(_TdiffphiT) , PdiffphiP(_PdiffphiP)
    { /* Constuctor for partial melt case. Constructor just stores important quantities*/ }
    Fluxes( T const & _Re_crit , T const & _A , T const & _B , T const & _l_nu, T const & _K_cond , T const & _nab_ad ) : 
       Re_crit(_Re_crit) , A_vis(_A) ,  B_conv(_B) , l_nu(_l_nu) , K_cond(_K_cond) , nab_ad(_nab_ad)  
    { /* Constructor for fully solid/molten or no settling case. Just stores important quantities*/ }
    Fluxes() = default;
    T LConv(T const& x)
    {
        u = velocity(x);
        T L_conv = B_conv * x * u ;
        return L_conv;
    }
    T velocity(T const& x)
    {
        T u_vis = A_vis*x;
        T u0 = Re_crit/l_nu*(-1 - pow( 1 + 4*u_vis*l_nu/Re_crit ,1.0/2))/ 2 ;
        return -u_vis*Re_crit/(u0*l_nu);
    }
    T LCond(T const& x)
    {        
        return K_cond * (x + nab_ad);
    }
    T LMix(T const& x , T const& DeltaS)
    {
        return B_conv/Cp * u * DeltaS * ((x+nab_ad)*TdiffphiT + PdiffphiP);
    }
    static T LGrav(T &C_grav , T &phi , T &DeltaS)
    {
        return C_grav * F_phi(phi) * DeltaS;
    }
    T u;
    private:
        T Re_crit , A_vis , B_conv , l_nu , K_cond , nab_ad , DeltaS , Cp , TdiffphiT , PdiffphiP , C_grav ; 
        static T F_phi(T const& phi){
            if(phi < 0.5){
                return std::max(pow(phi,3)/(1-phi)/1000,5.0/7*pow(phi,5.5)*(1-phi));
            }else{
                return std::min(2.0/9*phi*(1-phi) ,5.0/7*pow(phi,5.5)*(1-phi));
            }    
            
        }
};