import numpy as np
import scipy.integrate as integrate
import scipy.optimize as opt

# fitted to Spera (2011)
T0 = 3000
m = 0.6	
rho0 = 1.0/12.949 #A per ATOM
K0 = 13.20 #[GPa]	
K0_diffP = 8.238 
E0 = -20.5953
gamma0 = +0.1899
gamma0_prime = -1.940
a1 = 6*gamma0
a2 = -12*gamma0 +36*gamma0**2 -18*gamma0_prime
b_dim = 4
b = [+0.9821 ,+0.615, +1.31, -3.0 ,  -4.1]# eV/atom
k_boltz = 1.38064852e-23 #m2 kg s-2 K-1
m_u = 1.66053906660e-27 #kg
n = 5
CVkin = 1.5 * k_boltz
S0 = 3000 #?????? the mystery continues
Molar_mass = 100.385
JeV = 1.602176565e-19
def VinetP(rho):
    x = (rho0/rho)**(1.0/3)
    eta = 1.5*(K0_diffP -1 )
    return 3*K0*x**-2*(1-x)*np.exp(eta*(1-x))

def f_therm(T):
    return (T/T0)**m -1
def f_prime(T):
    return m*(T/T0)**m / T
def b_RT(rho):
    b_tot = 0*rho
    for n in range(b_dim):
        b_tot = b_tot + b[n]*(rho0/rho -1)**n 
    return JeV * b_tot #J per atom
def b_prime(rho):
    b_p = 0*rho
    for n in range(b_dim-1):
        b_p = b_p + b[n+1]*rho*(n+1) * (rho0/rho -1)**(n)
    return 1e30 *JeV * b_p #J m-3
def gamma0S_func(f, T0S):
    return 1.0/6 * T0**2/T0S**2 * (a1 +a2*f) *(2*f+1)
def T0S_func(f):
    return T0*np.sqrt( 1 +a1*f + 0.5*a2*f**2)
def Cv_func(rho,T):
    return (b_RT(rho)*f_prime(T) + CVkin)  # J K-1 atom-1
def RTpress_P(rho,T):
    P0 = VinetP(rho)
    P_E = -b_prime(rho)*(f_therm(T)-f_therm(T0))
    f = 0.5*((rho/rho0)**(2.0/3)-1)
    T0S = T0S_func(f)
    P_S = b_prime(rho)/(m-1) * ( T*(f_prime(T)-f_prime(T0S)) - T0*(f_prime(T0)-f_prime(T0S))) + gamma0S_func(f,T0S)*Cv_func(rho,T0S) * (T-T0)*rho*1e30    
    return P0 + (P_E + P_S)/1e9
def RTpress_S(rho,T):
    f = 0.5*((rho/rho0)**(2.0/3)-1)
    T0S = T0S_func(f)
    S = S0 + (b_RT(rho)/(m-1) *(f_prime(T) - f_prime(T0S)) + CVkin * np.log(T/T0S) )*5/(Molar_mass*m_u)
    return S 


def EoSthings(P,T):
    if P[0] < P[-1]:
        P = np.flip(P)

    rhocalc = np.zeros(len(P))
    deltacalc = np.zeros(len(P))
  
    Scalc = np.zeros(len(P))
    Cpcalc = np.zeros(len(P))
    Cpcalc2 = np.zeros(len(P))
    Cvcalc = np.zeros(len(P))


    rho_min = rho0*0.5
    rho_max = rho0*100
    if T == 6000:
        rho_min = rho0*0.5
    
    for i in range(len(P)):
        if T > 3000:
            rho_min = rho0/2
        def RTpress_P_fixT(rho):
            return RTpress_P(rho,T) - P[i]
        
        print(T , rho_min , rho_max , rhocalc[i-1] , P[i] , RTpress_P_fixT(rho_min) , RTpress_P_fixT(rho_max))
        
        rhocalc[i] = opt.brentq(RTpress_P_fixT , rho_min, rho_max)
        
        Scalc[i] = RTpress_S(rhocalc[i] , T)
        
        dT = T*1e-5
        def RTpress_P_fixdT(rho):
            return RTpress_P(rho,T+dT) - P[i]
       
        rhodrho = opt.brentq(RTpress_P_fixdT , rho_min, rho_max)
        SdS = RTpress_S(rhodrho, T+dT)
        
        deltacalc[i] = -(rhodrho-rhocalc[i])/dT * T/rhocalc[i]

        Cpcalc[i] = T * (SdS - Scalc[i])/dT

        dP = P[i]*1e-5
        def P_dP(rho):
            return RTpress_P(rho,T) - (P[i]+dP)

        rhodrho = opt.brentq(P_dP, rhocalc[i]*0.9, rhocalc[i]*1.01)
        beta_T = 1/rhocalc[i]*(rhodrho - rhocalc[i])/dP

        Cvcalc[i] = Cv_func(rhocalc[i],T)*5/(Molar_mass*m_u) 
        Cpcalc2[i] = Cvcalc[i]+ deltacalc[i]**2/T / (rhocalc[i]* Molar_mass *m_u*1e30/5 *beta_T/1e9)

        rho_max = rhocalc[i]*1.1
        if P[i] < 1e-9 or (P[i] < 1e-8 and T>5000):
            rho_max = rho_max
            

        if deltacalc[i] < 0:
            sys.exit()

    if P[0] > P[-1]:
        #not actually flipping atm because needs to be P=
        rhocalc = np.flip(rhocalc)
        deltacalc = np.flip(deltacalc)
        Scalc = np.flip(Scalc)
        Cpcalc = np.flip(Cpcalc)
        Cpcalc2 = np.flip(Cpcalc2)
        Cvcalc = np.flip(Cvcalc)

        P = np.flip(P)
       

    return rhocalc  *1e30*Molar_mass*m_u/5 , deltacalc / (rhocalc  *1e30*Molar_mass*m_u/5 ) , Scalc , Cpcalc #, Cvcalc 