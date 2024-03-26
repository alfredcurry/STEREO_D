import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
from matplotlib import lines
import math

'''Script for plotting snapshots of the internal evolution of a planet's properties. One can also produce a series of images to create an animation. Currently you have to manually decide what timesteps to use.'''

#plt.style.use('paper_style.mplstyl')

#Physical constants
M_E = 5.972e24 #kg
R_E = 6.371e6 #m
G_Newt =  6.67428e-11 

def Earth_to_km(R):
    return R*R_E/1000
def km_to_Earth(R):
    return R*1000/R_E
def add_km_axis(ax ):
    secax = ax.secondary_yaxis('right', functions=(km_to_Earth , Earth_to_km ))
    secax.set_ylabel("Radius [R$_{\oplus}$]")

#figure properties
figrows = 2

fig , axs = plt.subplots(figrows, 2 , sharex = 'col', figsize = (12 , 8))
fig .subplots_adjust( hspace= 0.1 , wspace= 0.2)

lstyls =  list(lines.lineStyles.keys())
lstyls = lstyls[:4]

scalx = "linear"
scaly = "linear"
paramlabels = [ "Temperature, $T$ [K]" ,'Melt Fraction, $\phi$' ,  r'Density, $\rho$ [kg/m$^3$]' , r"Viscosity, $\eta$ [Pa s]" ]

#Planet properties

M = 0.15

iron_melting = np.loadtxt("eos/LiquidIron/Iron_melting.csv", delimiter = ',')
axs[0,0].plot(iron_melting[:,0],iron_melting[:,1] , linestyle = 'dashed' , color = 'k')

run = 'Fitted_Booth'
MLtype = "NB"
Mdot = 0
folder = "And11_Litasov/Booth_new"
#colors = np.full(len(ts), 'r')
Ps = 5e8
Linit = 7e5

Mcore = 0 #1.95e24
Rcore = 0 #3.5e6
Lcore = 0
num = 1
core_frac = 0.3

tau = 0
Teq = 2320
Res = 800
ftime = 0.001


ts = [ 0 , 100.04 , 1001 , 10850 ,1.019e+10]
colors = pl.cm.jet_r(np.linspace(0,1,len(ts)))


i=0

for t in ts:
    #For animation
    #fig , axs = plt.subplots(2, 2, sharex = 'col' , figsize = (8, 6))
    #fig .subplots_adjust( wspace= 0.5)
    for i2 in range(figrows):

        axs[0,i2].set_xscale(scalx)
        axs[1,i2].set_xscale(scalx)

        axs[0,i2].set_yscale(scaly)
        axs[1,i2].set_yscale(scaly)
    
        axs[0,i2].set_ylabel(paramlabels[i2])
        axs[1,i2].set_ylabel(paramlabels[i2+2])

        axs[figrows-1,i2].set_xlabel("Pressure, $P$ [GPa]")

        
    if run == 'Fitted_Booth':
        txtfile = "Results/"+folder+"/results{7}MgSiO3_E_M_E{1:.3g}_Res_{0}_ftime_{4:.2g}corefrac0.3_cond_4_ML_{9}L{2:.2g}_Teq{6:.0f}tau{8:.3g}FITPs{3:.2g}t={5:.5g}.txt".format( Res , M , Linit , Ps , ftime ,t , Teq , num , tau , MLtype)
        scalefile = "Results/"+folder+"/scale{7}MgSiO3_E_M_E{1:.3g}_Res_{0}_ftime_{4:.2g}corefrac0.3_cond_4_ML_{9}L{2:.2g}_Teq{6:.0f}tau{8:.3g}FITPs{3:.2g}t={5:.5g}.txt".format( Res , M , Linit , Ps , ftime, t, Teq , num , tau , MLtype)
        otherfile = "Results/"+folder+"/other{7}MgSiO3_E_M_E{1:.3g}_Res_{0}_ftime_{4:.2g}corefrac0.3_cond_4_ML_{9}L{2:.2g}_Teq{6:.0f}tau{8:.3g}FITPs{3:.2g}t={5:.5g}.txt".format( Res , M , Linit , Ps , ftime, t, Teq , num , tau , MLtype )
        fluxfile =  "Results/"+folder+"/fluxes{7}MgSiO3_E_M_E{1:.3g}_Res_{0}_ftime_{4:.2g}corefrac0.3_cond_4_ML_{9}L{2:.2g}_Teq{6:.0f}tau{8:.3g}FITPs{3:.2g}t={5:.5g}.txt".format( Res , M , Linit , Ps , ftime, t, Teq , num , tau , MLtype )
        
    print (txtfile)
    val = np.loadtxt(txtfile)
    scale = np.loadtxt(scalefile)
    other = np.loadtxt(otherfile)
    
    rho = val[:,5]
    
    g = G_Newt*Mcore/(scale[2]*val[:,2])**2

    r = (np.append([Rcore] , scale[2]*val[:,2]))/ R_E
    m = np.append(0 , val[:,0]*scale[0]/M_E) 
    r1 = ((val[:,2]*scale[2] - np.diff(scale[2]*val[:,2], prepend = Rcore)/2))/R_E
    
    #pressure
    P =  val[:,1]*scale[1] / 1e9

    #Temperature
    mark = "None"
    
    T = val[:,3]*scale[3]
    axs[0,0].plot(P , T , color = colors[i] , linestyle = lstyls[i%len(lstyls)] , marker = mark)
    axs[0,0].set_ylim(400, 4100) 

    T_liq = val[rho<7000,10]
    T_sol = val[rho<7000,9]

    axs[0,0].plot(P[rho<7000] , T_liq , color = 'k' , linewidth=2)
    axs[0,0].plot(P[rho<7000] , T_sol , color = 'k' , linewidth=2)
    
    axs[1,0].set_xlim(0,70)
    axs[1,1].set_xlim(0,35)
    
    if t < 1e4:
        displaytime = "{0:.0f} yrs".format(t)
    else:
        displaytime = "{0:.3g}".format(t/10**np.floor(np.log10(t))) + r"$\times$" + "10$^{{{0}}}$ yrs".format(int(np.floor(np.log10(t))))  
    #density
    axs[1,0].plot(P, rho, color = colors[i] , linestyle = lstyls[i%len(lstyls)] ,  label = displaytime) #rho
    axs[1,0].set_ylim(3000,10000)

    #viscosity
    visc = other[rho<7000,2]*rho[rho<7000]
    
    #melt fraction
    phi = val[rho<7000,8]
    
    axs[0,1].plot(P[rho<7000], phi , color = colors[i] , linestyle = lstyls[i%len(lstyls)] , marker = mark) 
    axs[0,1].set_yscale("linear")
    axs[0,1].set_ylim(-0.05,1.05) 
    
    axs[1,1].plot(P[rho<7000], visc , linestyle = lstyls[i%len(lstyls)] ,  color = colors[i]) #)
    axs[1,1].set_ylim(1e-2,1e31) 
    axs[1,1].set_yscale("log")

    
    #For animation
    #fig.suptitle("Time = {0:.4g} yrs".format(t) , fontsize=16 , y = 0.94)

    #fig1name = './Figures/Rocky/Animation/ML_NB_Run_core_frac{2}M{1}t={0:.0f}.jpg'.format(t, M , core_frac)
    #fig.savefig(fig1name, bbox_inches = 'tight', dpi = 600)
    #plt.close()
    i=i+1

axs[0,0].text(0.05,0.95 ,"(a)",  horizontalalignment='center', verticalalignment='top' , fontsize = 14 , transform=axs[0,0].transAxes)
axs[1,0].text(0.05,0.95 ,"(b)",  horizontalalignment='center', verticalalignment='top' , fontsize = 14 , transform=axs[1,0].transAxes)
axs[0,1].text(0.95,0.95 ,"(c)",  horizontalalignment='center', verticalalignment='top' , fontsize = 14 , transform=axs[0,1].transAxes)
axs[1,1].text(0.95,0.95 ,"(d)",  horizontalalignment='center', verticalalignment='top' , fontsize = 14 , transform=axs[1,1].transAxes)


fig.legend( loc = 'center left' , bbox_to_anchor = [0.9 , 0.5])

plt.show()
pr = input("make a fig?")

if pr == 'yes' or pr == 'y':
    
    figname = './Figures/Paper1/Evolve_Mdot{2}M{0}Tss{1}_P.pdf'.format(M , Teq , Mdot)
    fig.savefig(figname, bbox_inches = 'tight')

    print ("saved: " , figname )
else: print ("not saved")