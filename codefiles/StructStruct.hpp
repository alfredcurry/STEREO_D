#ifndef _STRUCTURE_STRUCT__H
#define _STRUCTURE_STRUCT__H

#include <iostream>
#include <Eigen/Dense>

struct P_P{
    /**
     * @brief Structure for holding the various properties of the planet. Most are self explanatory other than 'y' which is an array of the scaled variables P,r,T,L, (in that order). P,T are defined at cell centres, r,L at cell edges. 'scale' contains scales of these variables. 'ym' contains the unscaled variables defined at the cells centres. 'm' is the scaled mass at cell edges and 'mm' at cell centres
     * 
     */

    Eigen::ArrayXd m , mm , dm , dmm, rho, scale, Cp , del_rho , nabla, nabAd, nabNonAd, T0 , P0 , H , E_mass_loss , g , l_mix , nu , L_conv , L_cond , L_grav , L_mix , k_cond , phi , rho_s , rho_m , Cp_s , Cp_m , Cp_avg , del_rho_s , del_rho_m , del_rho_avg , DeltaS , DeltaV , u_conv;
    Eigen::ArrayXXd y , dy , ym;
    double n , mu_m , Tsurf[7] , Psurf[7] , conds[3] , dt , M , opac , M_core , R_core , L_core , R_core0 , R0 , tau , Teq , Teff , a_crystal , visc_m , T_BL , P_BL , dT_BL[4] , dP_BL[4] , dviscdphi , dviscdT , dviscdP , z_BL , M_dot , M_total , z_edge , M_edge , z_min , P_min , g_surf , rho_fac , solid_core_mass[3] , M_from_day , M_dot_from_day;
    char Tbound , Pbound , nabtype , BL , ML , MassLossType = '0', RadioType= '0';
    Eigen::ArrayXd K , diffrhoP, diffrhoT , diffdel_rhoP, diffdel_rhoT , diffCpP , diffCpT , Bavg , supernabla , GM_R4 , diffGM_R4 , diffphiP , diffphiT ;
    Eigen::ArrayXXd diffnabla;
    int Jcore = 0;
    void resize(int J)
    {
        int I = scale.size();
        m.resize(J) ; mm.resize(J) ; dm.resize(J-1) ; dmm.resize(J-1); rho.resize(J); Cp.resize(J) ; Cp_avg.resize(J) ; del_rho.resize(J) ; del_rho_avg.resize(J) ; nabla.resize(J); nabAd.resize(J); nabNonAd.resize(J); T0.resize(J) ; P0.resize(J) ; H.resize(J-1) ; GM_R4.resize(J-1) ; diffGM_R4.resize(J-1) ; g.resize(J) ; l_mix.resize(J) ; nu.resize(J) ; L_conv.resize(J) ; L_cond.resize(J) ; L_grav.resize(J) ; L_mix.resize(J) ; k_cond.resize(J) ;  
        y.resize(J,I) ; dy.resize(J-1,I) ; ym.resize(J,I) ; K.resize(J) ; diffrhoP.resize(J); diffrhoT.resize(J) ; diffdel_rhoP.resize(J); diffdel_rhoT.resize(J); diffCpP.resize(J); diffCpT.resize(J) ; diffnabla.resize(J,I) ; supernabla.resize(J) ; Bavg.resize(J-1) ; 
        phi.resize(J) ; diffphiP.resize(J) ; diffphiT.resize(J) ; rho_s.resize(J) ; rho_m.resize(J) ; Cp_s.resize(J) ; Cp_m.resize(J) ; del_rho_s.resize(J) ; del_rho_m.resize(J) ; DeltaS.resize(J) ; DeltaV.resize(J) ; u_conv.resize(J) ; E_mass_loss.resize(J-1) ;
    }
    P_P(int I, int J):
        m(J) , mm(J) , dm(J-1) , dmm(J-1), rho(J), scale(I), Cp(J) , Cp_avg(J) , del_rho(J) , del_rho_avg(J) , nabla(J), nabAd(J), nabNonAd(J), T0(J) , P0(J) , H(J-1) , GM_R4(J-1) , diffGM_R4(J-1) , g(J) , l_mix(J) , nu(J) , L_conv(J) , L_cond(J) , L_grav(J) , L_mix(J) , k_cond(J) ,  
        y(J,I) , dy(J-1,I) , ym(J,I) , K(J) , diffrhoP(J), diffrhoT(J) , diffdel_rhoP(J), diffdel_rhoT(J), diffCpP(J), diffCpT(J) , diffnabla(J,I) , supernabla(J) , Bavg(J-1) , 
        phi(J) , diffphiP(J) , diffphiT(J) , rho_s(J) , rho_m(J) , Cp_s(J) , Cp_m(J)  , del_rho_s(J) , del_rho_m(J) , DeltaS(J) , DeltaV(J) , u_conv(J), E_mass_loss(J-1)
        {}
};

#endif