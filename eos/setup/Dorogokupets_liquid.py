import numpy as np
import scipy.integrate as integrate
import scipy.optimize as opt
#from numba import jit

#for gamma iron - Dorogojupets (2017)
K0 = 83.7 #GPa
K0_diffP =  5.97
K0_diff2P = 0
theta0 = 263
gamma0 = 2.003
gamma_inf = 0
e0 = 198e-6 #K-1
g = 0.884
beta = 1.168
a_s = 2.12

N_avo = 6.0221408e+23

rho0 = 1.0/7.957 #molcm-3 OR 1.0/49.026 #A-3 
R_gas =  8.31446261815324 #JK-1mol-1 
k_boltz = 1.38064852e-23 #m2 kg s-2 K-1
n = 1
Cvm = 3*n*R_gas
Molar_mass = 55.845 * n
T0 = 1811 #K 

def Debye(x):
    return x**3/(np.exp(x)-1)
def D_integ(x):
    return 3.0/x**3 * integrate.quad(Debye,0,x)[0]

def gamma_func(rho):
    return gamma_inf + (gamma0-gamma_inf)*(rho0/rho)**beta
def theta_func(rho, T):
    x = rho0/rho
    return T*x**-gamma_inf*np.exp((gamma0-gamma_inf)/beta*(1-x**beta))

def VinetP(rho):
    x = (rho0/rho)**(1.0/3)
    eta = 1.5*(K0_diffP -1 )
    return 3*K0*x**-2*(1-x)*np.exp(eta*(1-x))

def electronic(rho):
    x = rho0/rho
    return e0*x**g 
def P_electronic(rho,T):
    return 1.5*n*R_gas*electronic(rho)*g*rho* T**2 * 1e-3
    
def Dorog_P(rho, T):
    
    P0 = VinetP(rho)
    Pe = P_electronic(rho,T) -  P_electronic(rho,T0)
    gamma  = gamma_func(rho)
    theta = theta_func(rho , theta0)
    
    Ptherm = Cvm * gamma * rho * (theta/(np.exp(theta/T)-1) - theta/(np.exp(theta/T0)-1)) * 1e-3
    return P0 + Ptherm + Pe #

def Dorog_Cv(rho,T):
    theta = theta_func(rho ,theta0)
    Cv = Cvm * (theta/T)**2 * np.exp(theta/T)/(np.exp(theta/T)-1)**2 + a_s*R_gas
    Cve = Cvm*electronic(rho)*T
    return Cv + Cve


def EoSthings(P,T):
    if P[0] < P[-1]:
        P = np.flip(P)

    rhocalc = np.zeros(len(P))
    delta_rhocalc = np.zeros(len(P))
  
    Scalc = np.zeros(len(P))
    Cpcalc = np.zeros(len(P))
    rho_min = rho0*0.5
    rho_max = rho0*4
    if T > 1000:
        rho_min = rho0*0.6
    #if T > 4000:
        #rho_min = rho0*0.5
    
    for i in range(len(P)):
        
        def Dorog_P_fixT(rho):
            return Dorog_P(rho,T) - P[i]
        
        print(T , rho_min , rho_max , P[i] , Dorog_P_fixT(rho_min) , Dorog_P_fixT(rho_max))
        if T > 4000 and P[i] < 20:
            rhocalc[i] = 8000.0  
        elif T > 5000 and P[i] < 30:
            rhocalc[i] = 8000.0   
        else:
            rhocalc[i] = opt.brentq(Dorog_P_fixT , rho_min, rho_max)
        
            dT = T*1e-5
            def Dorog_P_fixdT(rho):
                return Dorog_P(rho,T+dT) - P[i]
       
            rhodrho = opt.brentq(Dorog_P_fixdT , rho_min, rho_max)
        
            delta_rhocalc[i] = -(rhodrho-rhocalc[i])/dT * T/(rhocalc[i]*rhocalc[i])

            drho = rhocalc[i]*1e-5
            K_T = rhocalc[i] * (Dorog_P(rhocalc[i]+drho,T)-P[i])/drho
            Cpcalc[i] = Dorog_Cv(rhocalc[i],T) + delta_rhocalc[i]*delta_rhocalc[i]*K_T*rhocalc[i]/T *1e3
            rho_max = rhocalc[i]
            #rho_min = 0.1*rhocalc[i]


    if P[0] > P[-1]:
        #not actually flipping atm because needs to be P=
        rhocalc = np.flip(rhocalc)
        delta_rhocalc = np.flip(delta_rhocalc)
        Scalc = np.flip(Scalc)
        Cpcalc = np.flip(Cpcalc)
        P = np.flip(P)
       

    return rhocalc *Molar_mass *1000, delta_rhocalc /(Molar_mass *1000 ), Scalc/Molar_mass , Cpcalc/Molar_mass *1000