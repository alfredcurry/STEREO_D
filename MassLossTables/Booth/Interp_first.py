import os
import numpy as np
from numpy.lib.function_base import meshgrid
import scipy.interpolate as interp
import sys

'''Script that takes Richard Booth's mass loss rates (in npz file) and puts them in a format to be used by my 2d hermite interpolator'''

a_index = int(sys.argv[1])

log_m = np.linspace(-2 , -0.0 , 20)

basefile = "./MassLossTables/Booth/"
#basefile = "./Boundary/No_redist_T_Rlimited_Tss2145"

data = np.load(basefile+"mass_loss.npz")
rho_fac = data['rho_fac']

a = data['a'][a_index]
Tss = 2320 * (a/0.0134) ** -0.5
print(a , Tss)
folder = "a{0:.4f}Tss{1:.1f}".format(a,Tss)
path = basefile+folder
if(os.path.exists(path) == False):
        os.mkdir(path)
        os.mkdir(path+"/data")
        os.mkdir(path+"/grid")
        os.mkdir(path+"/coeff")

log_m_max = -0.0

print(rho_fac)
mdot_g_data = data['mdot_g'][a_index,data['log_m']<log_m_max,:]
log_m_data = data['log_m'][data['log_m']<log_m_max]
print(np.shape(mdot_g_data))
for row in range(np.shape(mdot_g_data)[0]):
        nan_row = np.isnan(mdot_g_data[row,:])
        if nan_row.any():
                print("nan for log m " , log_m_data[row])
                if np.sum(nan_row) == 1:
                        nan_index = np.where(nan_row == True)
                        print("rho_fac" , rho_fac[nan_index])
                        non_nan_index =  np.where(nan_row == False)
                        print(rho_fac[non_nan_index] ,  mdot_g_data[row,non_nan_index][0])
                        clean_interpolate = interp.PchipInterpolator(rho_fac[non_nan_index] , mdot_g_data[row,non_nan_index][0])
                        print("replacing")
                        mdot_g_data[row,nan_index] = clean_interpolate(rho_fac[nan_index])


print(mdot_g_data)
logmdot_func = interp.RectBivariateSpline( log_m_data , data['rho_fac'] , np.log10(mdot_g_data) )

logmdot = logmdot_func(log_m , rho_fac , dx = 0, dy = 0)
logmdot_dlog_m = logmdot_func(log_m , rho_fac , dx = 1, dy = 0)
logmdot_drho_fac = logmdot_func(log_m , rho_fac ,  dx = 0, dy = 1)
logmdot_dlog_mdrho_fac = logmdot_func(log_m , rho_fac ,  dx = 1, dy = 1)


np.savetxt((basefile+folder+"/grid/log_m{0:.0f}.txt" ).format(len(log_m)) , log_m)
np.savetxt((basefile+folder+"/grid/rho_fac{0:.0f}.txt").format(len(rho_fac)), rho_fac)
np.savetxt((basefile+folder+"/data/logmdot{0:.0f}x{1:.0f}.txt").format(len(log_m),len(rho_fac)),logmdot)
np.savetxt((basefile+folder+"/data/logmdot_dlog_m{0:.0f}x{1:.0f}.txt").format(len(log_m),len(rho_fac)), logmdot_dlog_m )
np.savetxt((basefile+folder+"/data/logmdot_drho_fac{0:.0f}x{1:.0f}.txt").format(len(log_m),len(rho_fac)), logmdot_drho_fac )
np.savetxt((basefile+folder+"/data/logmdot_dlog_mdrho_fac{0:.0f}x{1:.0f}.txt").format(len(log_m),len(rho_fac)), logmdot_dlog_mdrho_fac )
print("saved")