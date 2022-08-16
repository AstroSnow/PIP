import pipreadmods as PIPpy
import matplotlib.pyplot as plt
import numpy as np


fname='../Data/'
#fname='../simdata/nleveltest_rad/Data/'
#fname='../simdata/nleveltest_rad_thick_st'

ds=PIPpy.pipread(fname,5)

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

fig, axs = plt.subplots(2, 2,dpi=300)
axs[0,0].plot(ds['xgrid'],ds['ion'],label='ion')
axs[0,0].plot(ds['xgrid'],ds['rec'],label='rec')
axs[0,0].plot(ds['xgrid'],ds['ion_rad'],label='ion_rad')
axs[0,0].plot(ds['xgrid'],ds['rec_rad'],label='rec_rad')
axs[0,0].set_yscale('log')
axs[0,0].set_xscale('log')
axs[0,0].legend()

axs[0,1].plot(ds['xgrid'],ds['vx_p'],label='vx_p')
axs[0,1].plot(ds['xgrid'],ds['vx_n'],label='vx_n')
axs[0,1].set_xscale('log')
axs[0,1].legend()

axs[1,0].plot(ds['xgrid'],ds['pr_p']/ds['ro_p']*5.0/6.0,label='T_p')
axs[1,0].plot(ds['xgrid'],ds['pr_n']/ds['ro_n']*5.0/3.0,label='T_n')
axs[1,0].legend()
axs[1,0].set_xscale('log')

axs[1,1].plot(ds['xgrid'],ds['nexcite1'],label='n0')
axs[1,1].plot(ds['xgrid'],ds['nexcite2'],label='n1')
axs[1,1].plot(ds['xgrid'],ds['nexcite3'],label='n2')
axs[1,1].plot(ds['xgrid'],ds['nexcite4'],label='n3')
axs[1,1].plot(ds['xgrid'],ds['nexcite5'],label='n4')
axs[1,1].plot(ds['xgrid'],ds['nexcite6'],label='c')
axs[1,1].legend()
axs[1,1].set_xscale('log')
axs[1,1].set_yscale('log')

plt.savefig('test_plot.png',dpi=300)

"""
fig, axs = plt.subplots(2, 2)
axs[0,0].plot(ds['xgrid'],ds['colrat'][:,1,6],label='colrat16')
axs[0,0].plot(ds['xgrid'],ds['colrat'][:,2,6],label='colrat26')
axs[0,0].plot(ds['xgrid'],ds['colrat'][:,3,6],label='colrat36')
axs[0,0].plot(ds['xgrid'],ds['colrat'][:,4,6],label='colrat46')
axs[0,0].plot(ds['xgrid'],ds['colrat'][:,5,6],label='colrat56')
axs[0,0].set_yscale('log')
axs[0,0].set_xscale('log')
axs[0,0].legend()

axs[1,0].plot(ds['xgrid'],ds['colrat'][:,6,1],label='colrat61')
axs[1,0].plot(ds['xgrid'],ds['colrat'][:,5,1],label='colrat51')
axs[1,0].plot(ds['xgrid'],ds['colrat'][:,4,1],label='colrat41')
axs[1,0].plot(ds['xgrid'],ds['colrat'][:,3,1],label='colrat31')
axs[1,0].plot(ds['xgrid'],ds['colrat'][:,2,1],label='colrat21')
axs[1,0].set_yscale('log')
axs[1,0].set_xscale('log')
axs[1,0].legend()
"""
