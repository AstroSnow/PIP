#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 13:58:39 2022

@author: ben
"""
import pipreadmods as PIPpy
import matplotlib.pyplot as plt
import numpy as np


#fname='../Data/'
#fname='../simdata/nleveltest_rad/Data/'
#fname='../simdata/nleveltest_rad_thick_st'
#fname='../simdata/T_6220_long/Data/' ;xmin=0.075;xmax=0.08
fname='../simdata/isca_midc_ltest/Data/' ;xmin=0.052;xmax=0.074

ds=PIPpy.pipread(fname,30)
dsm=PIPpy.pipread('../simdata/MHD_ref/',50)

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

fig, axs = plt.subplots(1, 4,dpi=300)
fig.set_size_inches(19.8, 4)
lthick=2.5

#axs[0,0].plot(ds['xgrid']/ds['time'],ds['vx_n']-ds['vx_p'],label='vx_n-vx_p',color='k',linewidth=lthick)
axs[0].plot(ds['xgrid']/ds['time'],ds['vx_n'],label='vx_n',color='r',linewidth=lthick)
axs[0].plot(ds['xgrid']/ds['time'],ds['vx_p'],label='vx_p',color='b',linewidth=lthick)
#axs[0,0].set_xscale('log')
axs[0].set_xlim([xmin,xmax])
#axs[0,0].set_ylim([-0.16,0.16])
axs[0].legend()
axs[0].set_xlabel('$x/t$')
axs[0].set_ylabel('Velocity')

axs[1].plot(ds['xgrid']/ds['time'],ds['pr_p']/ds['ro_p']*5.0/6.0,label='T_p',color='b',linewidth=lthick)
axs[1].plot(ds['xgrid']/ds['time'],ds['pr_n']/ds['ro_n']*5.0/3.0,label='T_n',color='r',linewidth=lthick)
axs[1].plot(dsm['xgrid']/dsm['time'],dsm['pr_p']/dsm['ro_p']*5.0/3.0,'k--',label='T (MHD)',linewidth=lthick)
#axs[0,1].set_yscale('log')
#axs[0,1].set_xscale('log')
axs[1].set_xlim([xmin,xmax])#axs[0,1].set_xlim([0.01,2.0])
axs[1].legend()
axs[1].set_xlabel('$x/t$')
axs[1].set_ylabel('Temperature')

axs[2].plot(ds['xgrid']/ds['time'],ds['ion'],label='ion',linewidth=lthick)
axs[2].plot(ds['xgrid']/ds['time'],ds['rec'],label='rec',linewidth=lthick)
axs[2].plot(ds['xgrid']/ds['time'],ds['ion_rad'],label='ion_rad',linewidth=lthick)
axs[2].plot(ds['xgrid']/ds['time'],ds['rec_rad'],label='rec_rad',linewidth=lthick)
axs[2].set_yscale('log')
#axs20].set_xscale('log')
axs[2].set_xlim([xmin,xmax])
axs[2].legend()
axs[2].set_xlabel('$x/t$')
axs[2].set_ylabel('$\Gamma$')

axs[3].plot(ds['xgrid']/ds['time'],ds['nexcite1'],label='n1',linewidth=lthick)
axs[3].plot(ds['xgrid']/ds['time'],ds['nexcite2'],label='n2',linewidth=lthick)
axs[3].plot(ds['xgrid']/ds['time'],ds['nexcite3'],label='n3',linewidth=lthick)
axs[3].plot(ds['xgrid']/ds['time'],ds['nexcite4'],label='n4',linewidth=lthick)
axs[3].plot(ds['xgrid']/ds['time'],ds['nexcite5'],label='n5',linewidth=lthick)
axs[3].plot(ds['xgrid']/ds['time'],ds['nexcite6'],label='c',linewidth=lthick)
axs[3].legend()
#axs[1,1].set_xscale('log')
axs[3].set_yscale('log')
axs[3].set_xlim([xmin,xmax])
axs[3].set_xlabel('$x/t$')
axs[3].set_ylabel('Density')

"""axs[1,1].plot(ds['xgrid']/ds['time'],ds['ro_p'],label='ro_p',color='b',linewidth=lthick)
axs[1,1].plot(ds['xgrid']/ds['time'],ds['ro_n'],label='ro_n',color='r',linewidth=lthick)
axs[1,1].plot(ds['xgrid']/ds['time'],ds['ro_n']+ds['ro_p'],label='ro_n+ro_p',color='g',linewidth=lthick)
axs[1,1].plot(dsm['xgrid']/dsm['time'],dsm['ro_p'],label='ro_p (MHD)',color='k',linewidth=lthick)
axs[1,1].legend()
axs[1,1].set_xscale('log')
#axs[1,1].set_yscale('log')
axs[1,1].set_xlim([0.01,2.0])
axs[1,1].set_xlabel('$x/t$')
axs[1,1].set_ylabel('Density')"""

plt.savefig('diracplot_shocks.png',dpi=300)
#plt.savefig('test_plot.png',dpi=300)

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
