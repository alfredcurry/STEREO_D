import numpy as np

#from RTpress import *
from Dorogokupets_liquid import *
basefile = "./eos/LiquidIron"

'''Script for tabulating equation of states as functions of P,T from EoS models. The eos must be specified above in imports, as well as the folder to save it into.'''


if basefile ==  "./eos/Stixrude":
    T = np.linspace(20 , 6000 , 300) #liquid: 1400 8000 #solid: 200 5000 LINEAR 1400/200 - 6000 #iron 20, 4100,
    Plow = np.linspace(0 , 122 , 300) #liquid: -10 1000 #solid -10 800 LINEAR 0 -300
    Phigh = np.linspace(123 , 900 , 600) #liquid: -10 1000 #solid -10 800 LINEAR 0 -300
    P = np.concatenate((Plow,Phigh))
else:
    T = np.linspace(1200 , 6000 , 300) #liquid: 1400 8000 #solid: 200 5000 LINEAR 1400/200 - 6000 #iron 20, 4100,
    P = np.linspace(10 , 400 , 1000)

Tgrid , Pgrid  =  np.meshgrid(T,P)

dT , dP = np.meshgrid(1e-3*T,1e-3*P)

rho = np.zeros((len(P),len(T)))
rho_dP = np.zeros((len(P),len(T)))
rho_dT = np.zeros((len(P),len(T)))
rho_dPdT = np.zeros((len(P),len(T))) 

delta = np.zeros((len(P),len(T)))
del_dP = np.zeros((len(P),len(T)))
del_dT = np.zeros((len(P),len(T)))
del_dPdT = np.zeros((len(P),len(T))) 

Cp = np.zeros((len(P),len(T)))
Cp_dP = np.zeros((len(P),len(T)))
Cp_dT = np.zeros((len(P),len(T)))
Cp_dPdT = np.zeros((len(P),len(T))) 

for i in range(len(T)) :
    if basefile == "./eos/Stixrude":
        EoS = Stixrude_EoS("Mg_pv")
        rho[:len(Plow),i] , delta[:len(Plow),i]  , Scalc , Cp[:len(Plow),i]  = EoS.EoSthings(Plow,T[i])
        EoS = Stixrude_EoS("Mg_ppv") #now including the pv-ppv phase transition 
        rho[len(Plow):,i] , delta[len(Plow):,i]  , Scalc , Cp[len(Plow):,i]  = EoS.EoSthings(Phigh,T[i])
    else:
        rho[:,i] , delta[:,i]  , Scalc , Cp[:,i]  = EoSthings(P,T[i])


np.savetxt((basefile+"/data/P{0:.0f}.txt").format(len(P)), P)
np.savetxt((basefile+"/data/T{0:.0f}.txt").format(len(T)) , T )

np.savetxt((basefile+"/data/rho{0:.0f}x{1:.0f}.txt").format(len(P),len(T)),rho)
np.savetxt((basefile+"/data/del_rho{0:.0f}x{1:.0f}.txt").format(len(P),len(T)),delta)
np.savetxt((basefile+"/data/Cp{0:.0f}x{1:.0f}.txt").format(len(P),len(T)),Cp)
