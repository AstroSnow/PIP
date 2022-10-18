#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 11:33:51 2022

@author: ben
"""
import pipreadmods as PIPpy
import matplotlib.pyplot as plt
import numpy as np


#fname='../Data/'
#fname='../simdata/nleveltest_rad/Data/'
#fname='../simdata/nleveltest_rad_thick_st'
fname='../simdata/T_6220_long/Data/' ;xmin=0.075;xmax=0.08; T0=6220.0
#fname='../simdata/isca_midc_ltest/Data/' ;xmin=0.052;xmax=0.074

ds=PIPpy.pipread(fname,30)
dsm=PIPpy.pipread('../simdata/MHD_ref/',50)

nntot=ds['nexcite1']+ds['nexcite2']+ds['nexcite3']+ds['nexcite4']+ds['nexcite5']

ndif=ds['ro_n']-nntot

#############################################################
#Temperature
T1=ds['pr_p'][-1]/ds['ro_p'][-1]*5.0/6.0
T1m=dsm['pr_p'][-1]/dsm['ro_p'][-1]*5.0/3.0

Tmhd=dsm['pr_p'][750]/dsm['ro_p'][750]*5.0/3.0*T0/T1m/1000.0

#############################################################
#Estimate the shock frame
b=np.argmin(np.gradient(ds['vx_p'])) #plasma shock location
vs=ds['xgrid'][b]/ds['time'] #plasma shock velocity
xs=ds['xgrid']/ds['time']-vs #shock grid
bn=np.argmin(np.gradient(ds['vx_n'])) #neutral shock location
vsn=ds['xgrid'][bn]/ds['time'] #neutral shock velocity

############################################################
#Frictional heating
fheat=0.5*ds['ac']*ds['ro_p']*ds['ro_n']*(np.abs(ds['vx_n']-ds['vx_p'])+np.abs(ds['vy_n']-ds['vy_p']))**2
#Thermal damping
tdamp=3.0/2.0*ds['ac']*ds['ro_p']*ds['ro_n']*(ds['pr_n']/ds['ro_n']-ds['pr_p']/ds['ro_p']/2.0)

fig, axs = plt.subplots(3, 1,dpi=300)
fig.set_size_inches(6, 9)
lthick=2.5

axs[0].plot(xs,fheat,label='T_p',color='k',linewidth=lthick)
axs[0].plot(xs,tdamp,color='r',linewidth=lthick)
axs[0].plot(xs,-tdamp,color='b',linewidth=lthick)
#axs[0].plot(xs,ds['pr_n']/ds['ro_n']*5.0/3.0*T0/T1/1000.0,label='T_n',color='r',linewidth=lthick)
#axs[0].plot([-1,1],[Tmhd,Tmhd],'k--',label='T (MHD)',linewidth=lthick)
axs[0].set_yscale('log')
#axs[0,1].set_xscale('log')
axs[0].set_xlim([xmin-vs,xmax-vs])#
axs[0].set_ylim([1.0e-5,1.0e-1])
#axs[0].legend()
axs[0].set_xlabel('$x_s$')
axs[0].set_ylabel('Power')

axs[1].plot(xs,ds['pr_p']/ds['ro_p']*5.0/6.0*T0/T1/1000.0,label='T_p',color='b',linewidth=lthick)
axs[1].plot(xs,ds['pr_n']/ds['ro_n']*5.0/3.0*T0/T1/1000.0,label='T_n',color='r',linewidth=lthick)
axs[1].plot([-1,1],[Tmhd,Tmhd],'k--',label='T (MHD)',linewidth=lthick)
#axs[0,1].set_yscale('log')
#axs[0,1].set_xscale('log')
axs[1].set_xlim([xmin-vs,xmax-vs])#axs[0,1].set_xlim([0.01,2.0])
axs[1].legend()
axs[1].set_xlabel('$x_s$')
axs[1].set_ylabel('Temperature [kK]')

axs[2].plot(xs,ds['ion_loss'],label='P_n',color='k',linewidth=lthick)
#axs[1,1].plot(xs,ds['pr_p'],label='P_p',color='b',linewidth=lthick)
#axs[1,1].plot(xs,ds['pr_p']+ds['pr_n'],label='P_p+P_n',color='g',linewidth=lthick)
#axs[0,0].set_xscale('log')
axs[2].set_xlim([xmin-vs,xmax-vs])
#axs[0,0].set_ylim([-0.16,0.16])
#axs[1,1].legend()
axs[2].set_xlabel('$x_s$')
axs[2].set_ylabel('loss rate')

plt.savefig('heatterms_plot_upc.png',dpi=300)


