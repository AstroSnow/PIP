import pipreadmods as PIPpy
import matplotlib.pyplot as plt
import numpy as np


fname='../Data/'
#fname='../simdata/nleveltest_rad/Data/'
#fname='../simdata/nleveltest_rad_thick_st'
#fname='../simdata/col_t10000/' ;sname='col_10000' #xi_n=0.94787198761607150
#fname='../simdata/rad_thin_t10000/' ;sname='rad_thin_10000' #xi_n=0.94787198863788269
#fname='../simdata/rad_thick_t10000/' ;sname='rad_thick_10000' #xi_n=0.94787198863788269

ds=PIPpy.pipread(fname,0)

#ds2=PIPpy.pipread('../simdata/rad_thick_t10000/',10)

nntot=ds['nexcite1']+ds['nexcite2']+ds['nexcite3']+ds['nexcite4']+ds['nexcite5']

ndif=ds['ro_n']-nntot

#plt.plot(ds['pr_p']/ds['ro_p'])
#plt.plot(ds['vx_n'])
#plt.plot(ds['by'])
#plt.plot(ds['nexcite4']/nntot)
#plt.plot(ds['ion']/ds['ion'][-1])
#plt.plot(ds['rec']/ds['rec'][-1])
#plt.plot(-ds['rec']*ds['ro_p']+ds['ion']*ds['ro_n'])
#plt.plot(nntot)
#plt.plot(ds['ion_loss'])
#plt.plot(ds['ion_rad']/ds['ion_rad'][-1])
#plt.plot(ds['rec_rad']/ds['rec'])

fig, axs = plt.subplots(2, 1)
axs[0].plot(ds['xgrid'],ds['vy_n'],label='ion')
axs[0].plot(ds['xgrid'],ds['vy_p'],label='rec')
axs[0].plot(ds['xgrid'],ds['vz_n'],label='ion_rad')
axs[0].plot(ds['xgrid'],ds['vz_p'],label='rec_rad')
axs[0].set_yscale('log')
#axs[0,0].set_xscale('log')
axs[0].legend()

axs[1].plot(ds['xgrid'],ds['vx_n'],label='vx_n')
axs[1].set_yscale('log')
axs[1].legend()