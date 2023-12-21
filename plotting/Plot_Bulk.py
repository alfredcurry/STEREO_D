import numpy as np
from numpy.lib.function_base import meshgrid
import scipy.interpolate as interp
import matplotlib.pyplot as plt
import matplotlib.ticker

'''Script for plotting the evolution of some bulk properties (including crystalisation) for a planet.'''

#plt.style.use('paper_style.mplstyl')

M_E = 5.972e24 #kg
R_E = 6.371e6 #m
G_Newt =  6.67428e-11 

M_dot = 0
Mcore = 0

Lcore = 0
tau = 0
ftime = 0.001
M = 0.15
thisM=M
fig_version = "all"
Pss = [5e8]
num = 1
Linit = 7e5
Res = 800

Tmin = 0
Tmax = 4000
tmax = None

fig_therm , ax_therm  = plt.subplots(3,1 , figsize = (6, 10) , sharex = True , gridspec_kw={'height_ratios': [1, 1.2, 1.5]})#

i=0
mark = [ 'None', 'None' , 'None', 'None', 'None', 'None']
lstyl = ['solid' , 'dashed' , 'dotted' ]
colours = ['b' , 'k' , 'y' , 'g' , 'y' , 'cyan']

MLtype = "NB"
Teq = 2320

folder = "And11_Litasov/Booth_new"
logLsize = 300
gsize = 40
i = 0
for Ps in Pss:
    if M_dot == 1:
        bulkfile = "Results/"+folder+"/bulk{6}MgSiO3_E_M_E{1:.3g}_Res_{0}_ftime_{4:.2g}corefrac0.3_Booth_M_dot_cond_4_ML_{8}L{2:.2g}_Teq{5:.0f}tau{7:.3g}FITPs{3:.2g}.txt".format( Res , M , Linit , Ps , ftime , Teq , num , tau , MLtype)
    elif M_dot == 0: 
        bulkfile = "Results/"+folder+"/bulk{6}MgSiO3_E_M_E{1:.3g}_Res_{0}_ftime_{4:.2g}corefrac0.3_cond_4_ML_{8}L{2:.2g}_Teq{5:.0f}tau{7:.3g}FITPs{3:.2g}.txt".format( Res , M , Linit , Ps , ftime , Teq , num , tau , MLtype)
    else:
        bulkfile = "Results/"+folder+"/bulk{8}MgSiO3_E_M_E{1:.3g}_Res_{0}_ftime_{4:.2g}corefrac0.3_PBC_M_dot_cond_4L{2:.2g}_Teq{5:.0f}tau{6:.3g}TABPs{3:.2g}.txt".format( Res , M , Linit , Ps , ftime , Teq , tau ,  M_dot , num)
    
    print(bulkfile)
    bulk = np.loadtxt(bulkfile)
    m_row = 14+3
    bulk =bulk[:-1,:]
    m = bulk[bulk[:,m_row] > 0.1*bulk[0,m_row], m_row]
    last_index = len(m)
    print(last_index)
    #m = bulk[:last_index,m_row+1] 
    m_tot = bulk[:last_index,m_row+1]
    #print(m_tot)
    if M == 0.03:
        m_core = bulk[:last_index,m_row+4] 
    elif M == 0.1:
        m_core = np.full(last_index , 1.95e+23)
    r = bulk[:last_index,2]
    z_edge = bulk[:last_index,m_row+3] 

    r_tot = r+z_edge
    Rcore = bulk[:last_index,m_row+5]

    m_tot = bulk[:last_index,m_row+1]
    m_dot = bulk[:last_index,m_row+2]/M_E*(365*24*3600*1e9)
    m_core = bulk[:last_index,m_row+4] 

    g = G_Newt*m_tot/(r+z_edge)**2
    print("g" , g)
    L = bulk[:last_index,4]
    Lcore = bulk[:last_index,m_row+6]
    L_tild = L/(2*np.pi*r**2)
    t = bulk[:last_index,0]

    #Tn = Tnfunc(np.log10(L_tild))
    
    Tc = bulk[:last_index,3]
    Ts = bulk[:last_index,6]
    T_CMB = bulk[:last_index,m_row+7]
    
  #  ax_therm[0].plot(t, m[:]/ M_E, color = colours[i] , linestyle = "dashed")
 #   ax_therm[0].plot(t, m_tot[:]/ M_E, color = colours[i] , linestyle = "solid")
#    ax_therm[0].plot(t, m_core[:]/ M_E, color = colours[i] , linestyle = "dotted")
    
    ax_therm[0].plot(t , L/(4*np.pi*r_tot**2) , color = colours[i] , label = "$P_0$ = {0} GPa".format(Ps/1e9) , linestyle = "solid" )
    #ax_therm[0].plot(t , Lcore , color = colours[i] , linestyle = "solid" )


    ax_therm[1].plot(t , Tc , color = colours[i] , marker = mark[i] , linestyle = "dotted" , label = "Central Temperature")
    ax_therm[1].plot(t , T_CMB , color = colours[i] , marker = mark[i] , linestyle = "dashed" , label = "Temperature at core-mantle boundary")
    ax_therm[1].plot(t , Ts , color = colours[i] , marker = mark[i] , linestyle = "solid" , label = "Temperature at 1 GPa")
    
    #ax_therm[0].plot(t, m[:]/ M_E, color = colours[i] , linestyle = "dashed")
    #ax_therm[0].plot(t, m_tot[:]/ M_E, color = colours[i] , linestyle = "solid")
    #ax_therm[0].plot(t, m_core[:]/ M_E, color = colours[i] , linestyle = "dotted")
    
    def R_plot(ax):
        ax.plot(t, r_tot/R_E , color = 'k' , label = 'Total radius')
        ax.plot(t, r/R_E , color = colours[i] , label = 'Radius at 1 GPa' )
        ax.plot(t, Rcore/R_E , color = 'r'  , label = 'Core radius')

        #ax.fill_between(t, r/R_E , Rcore/R_E , color = 'r' , alpha = 0.5 )
        ax.plot(t, np.clip(bulk[:last_index,9]/R_E, Rcore/R_E, 1) ,linestyle = 'dashed', color = 'r' , label = 'Liquidus', zorder=0)
        #ax.fill_between(t, np.clip(bulk[:last_index,9]/R_E, Rcore/R_E, 1) ,  Rcore/R_E , color = 'b' , alpha = 0.2 )
        #ax.fill_between(t, np.clip(bulk[:last_index,10]/R_E, Rcore/R_E, 1) ,   Rcore/R_E , color = 'b' , alpha = 0.2 )
        #ax.fill_between(t, np.clip(bulk[:last_index,8]/R_E, Rcore/R_E, 1) ,   Rcore/R_E , color = 'b' , alpha = 0.3 )
        ax.plot(t, np.clip(bulk[:last_index,10]/R_E, Rcore/R_E, 1) ,linestyle = 'dashdot' , color = colours[i] , label = '$\phi = \phi_c$', zorder=0)

        r_small_up = bulk[bulk[:,13]<=bulk[:,12],12]
        t_small = bulk[bulk[:,13]<=bulk[:,12],0]
        r_small_down = bulk[bulk[:,13]<=bulk[:,12],13]
        r_c_small = bulk[bulk[:,13]<=bulk[:,12],m_row+5]
            
        r_small_down1 = r_small_down[(r_small_down>r_c_small) | (t_small <1e6)]
        r_small_up = r_small_up[(r_small_down>r_c_small) | (t_small <1e6)]
        t_small = t_small[(r_small_down>r_c_small) | (t_small <1e6)]
            
        r_small_down = r_small_down1
         
        ax.plot(t_small[:], r_small_down[:]/R_E , linestyle =  "dashed", color = 'b' , zorder=0)
        ax.plot(t_small[:], r_small_up[:]/R_E , linestyle =  "dashed", color = 'b' , label = '$\phi = $0.05' , zorder=0 )
        ax.plot(t[bulk[:,11]<bulk[:,8]], bulk[bulk[:,11]<bulk[:,8],8]/R_E , linestyle =  "dotted", color = colours[i] , label = 'Solidus', zorder=0)
        ax.plot(t[bulk[:,11]<bulk[:,8]], bulk[bulk[:,11]<bulk[:,8],11]/R_E ,linestyle =  "dotted", color = colours[i] , zorder=0)
        t_solid = 1e10
        #ax.fill_between(t , r_tot  , r/R_E , color = 'r' , alpha = 0.7 )
        #ax.fill_between(t[t>=t_solid] , (r_tot[t>=t_solid] - 10e3)/R_E , r[t>=t_solid]/R_E , color = 'b' , alpha = 0.7 )
        #ax.fill_between(t[t>=t_solid] , r_tot[t>=t_solid]/R_E , (r_tot[t>=t_solid] - 10e3)/R_E  , color = 'r' , alpha = 0.7 )


    
    R_plot(ax_therm[2])
 
    i = i+1

ax_therm[0].set_ylabel("Mean flux [Wm$^{-2}$]")
ax_therm[0].set_yscale("log")
ax_therm[0].set_yticks(np.logspace(-3,6,10))
#ax_therm[1].set_ylim(9e11,2e20)
#ax_therm[0].set_ylabel("Mass [M$_\oplus$]")

ax_therm[1].set_ylabel("Temperature [K]")
ax_therm[1].set_ylim(Tmin,Tmax)

ax_therm[len(ax_therm)-1].set_xlabel("Time [yrs]")
ax_therm[len(ax_therm)-1].set_xscale("log")
ax_therm[len(ax_therm)-1].set_xlim(100 , tmax)


locmaj = matplotlib.ticker.LogLocator(base=10,numticks=12) 
locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=12)

ax_therm[len(ax_therm)-1].xaxis.set_major_locator(locmaj)
ax_therm[len(ax_therm)-1].xaxis.set_minor_locator(locmin)
ax_therm[len(ax_therm)-1].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

#ax_therm[len(ax_therm)-1].xaxis.grid(True, which='minor')
def Earth_to_km(R):
    return R*R_E/1000
def km_to_Earth(R):
    return R*R_E/1000

def SI_to_CGS_flux(F):
    return F/1e-7*1e-4
def CGS_to_SI_flux(F):
    return F*1e-7/1e-4

def config_R_plot(ax):
    #ax.set_title("Radius [R$_\oplus$]" , loc='left')
    ax.set_ylabel("Radius [R$_\oplus$]" )
    ax.set_xlabel("Time [yrs]")
    secax = ax.secondary_yaxis('right', functions=(Earth_to_km , km_to_Earth))
    #ax.set_title("Radius [km]", loc='right')
    secax.set_ylabel("Radius [km]" )
    ax.xaxis.set_major_locator(locmaj)
    ax.xaxis.set_minor_locator(locmin)
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())


config_R_plot(ax_therm[2])

#dummy_lines = []
#lnstyl = ["solid" , "dashed" , "dotted" ]
#typelegend = ["Total" , "1 GPa" , "Core"]

#for j in range(len(typelegend)):
 #   dummy_lines.append(ax_mass[0].plot([],[], c="black", ls = lnstyl[j])[0])
#legend2 = ax_mass[0].legend([dummy_lines[i] for i in range(len(lnstyl))], typelegend,loc = 3 , bbox_to_anchor=(1, 0.5))
#ax_mass[0].add_artist(legend2, )

#ax_therm[0].legend()
secLax = ax_therm[0].secondary_yaxis('right', functions=(SI_to_CGS_flux , CGS_to_SI_flux))
secLax.set_ylabel("Mean flux [erg s$^{-1}$cm$^{-2}$]")
#secLax.set_yticks(np.logspace(0,9,6))
#secLax.set_yticklabels(SI_to_CGS_flux(np.logspace(-3,6,5)))

ax_therm[1].legend(loc = "lower left" , bbox_to_anchor=(0, 0))
ax_therm[2].legend(loc = "lower right" , bbox_to_anchor=(1, 0.05))

#plt.close(fig_therm)

plt.show()

pr = input('make a fig?')
if pr == 'yes' or pr == 'y':
    
    figname = './Figures/Paper1/Bulk_Mdot{2}_M{0}Tss{1}.pdf'.format(thisM , Teq , M_dot)

    fig_therm.savefig(figname, bbox_inches = 'tight', dpi = 1000)
    
    print ("saved: " , figname )
else: print ("not saved")