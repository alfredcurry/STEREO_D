import numpy as np
from numpy.lib.function_base import meshgrid
import scipy.interpolate as interp
import matplotlib.pyplot as plt

''' Script for taking a tabulated grid on P,T for an EoS (produced by Setup_EoS) and smoothly interpolating it onto a second grid, with derivatives. This is useful because the 2D hermite interpolation requires derivatives'''

fig , ax = plt.subplots(2,2, sharex=True , sharey=True)
fig2 , ax2 = plt.subplots(2,2, sharex=True , sharey=True)

T = np.linspace( 1000 ,8000 , 300) 
P = np.linspace(20 , 3000 , 800)

basefile = "./eos/EpsilonIron"
Psize = 1000
Tsize = 300

rhogrid = np.loadtxt((basefile+"/data/rho{0}x{1}.txt").format(Psize , Tsize))
Cpgrid = np.loadtxt((basefile+"/data/Cp{0}x{1}.txt").format(Psize , Tsize))
del_rhogrid = np.loadtxt((basefile+"/data/del_rho{0}x{1}.txt").format(Psize , Tsize))
Pgrid = np.loadtxt((basefile+"/data/P{0:.0f}.txt").format(Psize))
Tgrid = np.loadtxt((basefile+"/data/T{0:.0f}.txt").format(Tsize))

rhofunc = interp.RectBivariateSpline(Pgrid ,Tgrid , rhogrid )
Cpfunc = interp.RectBivariateSpline(Pgrid ,Tgrid , Cpgrid )
del_rhofunc = interp.RectBivariateSpline( Pgrid ,Tgrid , del_rhogrid )

rho = rhofunc(P , T , dx = 0, dy = 0)
rho_dP = rhofunc(P , T , dx = 1, dy = 0)
rho_dT = rhofunc(P ,T ,  dx = 0, dy = 1)
rho_dPdT = rhofunc(P , T ,  dx = 1, dy = 1)

Cp = Cpfunc(P ,T ,  dx = 0, dy = 0)
Cp_dP = Cpfunc(P ,T ,  dx = 1, dy = 0)
Cp_dT = Cpfunc(P ,T ,  dx = 0, dy = 1)
Cp_dPdT = Cpfunc(P , T ,  dx = 1, dy = 1)

del_rho = del_rhofunc(P ,T ,  dx = 0, dy = 0)
del_rho_dP = del_rhofunc(P , T , dx = 1, dy = 0)
del_rho_dT = del_rhofunc(P , T ,  dx = 0, dy = 1)
del_rho_dPdT = del_rhofunc( P , T , dx = 1, dy = 1)

print (np.shape(rho))
print (np.shape(rhogrid))

c1 = ax[0,0].contourf( T , P ,  rho) 
c2 = ax[0,1].contourf(T , P , Cp) 
c3 = ax[1,0].contourf( T , P , del_rho)

d1 = ax2[0,0].contourf(Tgrid , Pgrid , rhogrid) 
d2 = ax2[0,1].contourf( Tgrid , Pgrid , Cpgrid) 
d3 = ax2[1,0].contourf( Tgrid , Pgrid , del_rhogrid ) 
fig.colorbar(c1 , ax = ax[0,0])
fig.colorbar(c2 , ax = ax[0,1])
fig.colorbar(c3, ax = ax[1,0])
fig2.colorbar(d1 , ax = ax2[0,0])
fig2.colorbar(d2 , ax = ax2[0,1])
fig2.colorbar(d3, ax = ax2[1,0])

plt.show()

np.savetxt((basefile+"/grid/P{0:.0f}.txt" ).format(len(P)) , P)
np.savetxt((basefile+"/grid/T{0:.0f}.txt").format(len(T)), T)
np.savetxt((basefile+"/data/{2:.0f}rho{0:.0f}x{1:.0f}.txt").format(len(P),len(T),Psize),rho)
np.savetxt((basefile+"/data/{2:.0f}rho_dP{0:.0f}x{1:.0f}.txt").format(len(P),len(T),Psize), rho_dP )
np.savetxt((basefile+"/data/{2:.0f}rho_dT{0:.0f}x{1:.0f}.txt").format(len(P),len(T),Psize), rho_dT )
np.savetxt((basefile+"/data/{2:.0f}rho_dPdT{0:.0f}x{1:.0f}.txt").format(len(P),len(T),Psize), rho_dPdT )
np.savetxt((basefile+"/data/{2:.0f}del_rho{0:.0f}x{1:.0f}.txt").format(len(P),len(T),Psize),del_rho)
np.savetxt((basefile+"/data/{2:.0f}del_rho_dP{0:.0f}x{1:.0f}.txt").format(len(P),len(T),Psize), del_rho_dP )
np.savetxt((basefile+"/data/{2:.0f}del_rho_dT{0:.0f}x{1:.0f}.txt").format(len(P),len(T),Psize), del_rho_dT )
np.savetxt((basefile+"/data/{2:.0f}del_rho_dPdT{0:.0f}x{1:.0f}.txt").format(len(P),len(T),Psize), del_rho_dPdT )
np.savetxt((basefile+"/data/{2:.0f}Cp{0:.0f}x{1:.0f}.txt").format(len(P),len(T),Psize),Cp)
np.savetxt((basefile+"/data/{2:.0f}Cp_dP{0:.0f}x{1:.0f}.txt").format(len(P),len(T),Psize), Cp_dP )
np.savetxt((basefile+"/data/{2:.0f}Cp_dT{0:.0f}x{1:.0f}.txt").format(len(P),len(T),Psize), Cp_dT )
np.savetxt((basefile+"/data/{2:.0f}Cp_dPdT{0:.0f}x{1:.0f}.txt").format(len(P),len(T),Psize), Cp_dPdT )
