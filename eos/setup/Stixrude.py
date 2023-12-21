import numpy as np
import scipy.integrate as integrate
import scipy.optimize as opt
#from numba import jit

class Stixrude_EoS:
    def __init__(self , material):
        self.R_gas =  8.31446261815324 #JK-1mol-1 
        self.k_boltz = 1.38064852e-23 #m2 kg s-2 K-1
        
        if material == "Mg_pv":
            #for Mg perovskite - Stix-II
            self.K0 = 251 #247 #247 #
            self.K0_diffP =  4.1 #3.97#
            self.K0_diff2P = -0.016
            self.theta0 = 905
            self.gamma0 = 1.57
            self.q0 = 1.1
            self.a1 = 6*self.gamma0   
            self.a2 = -12*self.gamma0 +36*self.gamma0**2 -18*self.q0*self.gamma0
            self.rho0 = 1.0/24.45 #*(4.10/4.1057295971) #mol cm-3 1.0/168.27 #A-3 
            self.n = 5
            self.Cvm = 3*self.n*self.R_gas
            self.Molar_mass = 100.385
            self.T0 = 1873 #K #Sinogeikin (2004)

        elif material == "Mg_ppv":
            #for Mg post-perovskite - Stix-II (deep mantle)

            self.K0 = 231
            self.K0_diffP =  4.0
            self.theta0 = 855
            self.gamma0 = 1.89
            self.q0 = 1.1
            self.a1 = 6*self.gamma0   
            self.a2 = -12*self.gamma0 +36*self.gamma0**2 -18*self.q0*self.gamma0
            self.rho0 = 1.0/24.42
            self.n = 5
            self.Cvm = 3*self.n*self.R_gas
            self.Molar_mass = 100.385
            self.T0 = 2535 #K #Guignot (2007)
            
            
    def Debye(self , x):
        return x**3/(np.exp(x)-1)
    def D_integ(self , x):
        return 3.0/x**3 * integrate.quad(self.Debye,0,x)[0]

    def gamma_func(self , f , theta):
        return 1.0/6 * self.theta0**2/theta**2 * (self.a1 +self.a2*f) *(2*f+1)

    def theta_func(self , f, T):
        return T*np.sqrt(np.abs( 1 +self.a1*f + 0.5*self.a2*f**2))

    def Stix_P(self , rho, T):
    
        b1 = 9*self.K0
        b2 = 27*self.K0*(self.K0_diffP-4)
        b3 = 0 # 9*K0*(9*K0_diffP**2 - 63*K0_diffP + 9 *K0*K0_diff2P+143)
        f = 0.5*((rho/self.rho0)**(2.0/3)-1)
        theta = self.theta_func(f , self.theta0)
        gamma  = self.gamma_func(f , theta)
    
        Ptherm = self.Cvm * gamma * rho * (T*self.D_integ(theta/T) - self.T0*self.D_integ(theta/self.T0)) * 1e-3
        return 1.0/3 * (1+2*f)**2.5 * ( b1*f + 0.5*b2*f**2 + 1.0/6 * b3*f**3) + Ptherm

    def Stix_S(self , rho,T):
        f = 0.5*((rho/self.rho0)**(2.0/3)-1)
        theta = self.theta_func(f, self.theta0)
        return -self.Cvm *( np.log(1-np.exp(-theta/T)) - 4.0/3*self.D_integ(theta/T))


    def EoSthings(self,P,T):
        if P[0] < P[-1]:
            P = np.flip(P)

        rhocalc = np.zeros(len(P))
        deltacalc = np.zeros(len(P))
  
        Scalc = np.zeros(len(P))
        Cpcalc = np.zeros(len(P))
        rho_min = self.rho0*0.7
        rho_max = self.rho0*4
        if T >= 5000:
            rho_min = self.rho0*0.5
    
        for i in range(len(P)):
        
            def Stix_P_fixT(rho):
                return self.Stix_P(rho,T) - P[i]
        
            print(T , rho_min , rho_max , P[i] , Stix_P_fixT(rho_min) , Stix_P_fixT(rho_max))
        
            rhocalc[i] = opt.brentq(Stix_P_fixT , rho_min, rho_max)
        
            Scalc[i] = self.Stix_S(rhocalc[i] , T)
        
            dT = T*1e-5
            def Stix_P_fixdT(rho):
                return self.Stix_P(rho,T+dT) - P[i]
       
            rhodrho = opt.brentq(Stix_P_fixdT , rho_min, rho_max)
            SdS = self.Stix_S(rhodrho, T+dT)
        
            deltacalc[i] = -(rhodrho-rhocalc[i])/dT * T/(rhocalc[i]*rhocalc[i])

            Cpcalc[i] = T * (SdS - Scalc[i])/dT
        #rho_max = rhocalc[i]

        if P[0] > P[-1]:
            #not actually flipping atm because needs to be P=
            rhocalc = np.flip(rhocalc)
            deltacalc = np.flip(deltacalc)
            Scalc = np.flip(Scalc)
            Cpcalc = np.flip(Cpcalc)
            P = np.flip(P)
       

        return rhocalc * self.Molar_mass *1000, deltacalc /(self.Molar_mass *1000 ), Scalc/self.Molar_mass , Cpcalc/self.Molar_mass *1000
