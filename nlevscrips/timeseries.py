#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 15:23:28 2022

@author: ben
"""
import pipreadmods as PIPpy
import matplotlib.pyplot as plt
import numpy as np


#fname='../Data/'
#fname='../simdata/nleveltest_rad/Data/'
#fname='../simdata/nleveltest_rad_thick_st'
fname='../simdata/upc_time10/' ;xmin=0.075;xmax=0.08; T0=6220.0
#fname='../simdata/isca_midc_ltest/Data/' ;xmin=0.052;xmax=0.074

ds1=PIPpy.pipread(fname,1)
ds10=PIPpy.pipread(fname,10)
ds100=PIPpy.pipread('../simdata/upc_time100/',10)
dsm=PIPpy.pipread('../simdata/MHD_ref/',50)

fig, axs = plt.subplots(3, 4,dpi=300)
fig.set_size_inches(15, 8)
#fig.tight_layout(pad=2.0)
plt.subplots_adjust(wspace=0.3)
lthick=2.5

#axs[0,0].plot(ds['xgrid']/ds['time'],ds['vx_n']-ds['vx_p'],label='vx_n-vx_p',color='k',linewidth=lthick)
axs[0,0].plot(ds1['xgrid']/ds1['time'],ds1['vx_n'],label='vx_n',color='r',linewidth=lthick)
axs[0,0].plot(ds1['xgrid']/ds1['time'],ds1['vx_p'],label='vx_p',color='b',linewidth=lthick)
axs[0,0].set_xscale('log')
axs[0,0].set_xlim([0.001,4])
#axs[0,0].set_ylim([-0.16,0.16])
#axs[0,0].legend()
#axs[0,0].set_xlabel('$x/t$')
axs[0,0].set_ylabel('$v_x$')

T1=ds1['pr_p'][-1]/ds1['ro_p'][-1]*5.0/6.0
T1m=dsm['pr_p'][-1]/dsm['ro_p'][-1]*5.0/3.0
Tmhd=dsm['pr_p'][750]/dsm['ro_p'][750]*5.0/3.0*T0/T1m/1000.0
axs[0,1].plot(ds1['xgrid']/ds1['time'],ds1['pr_p']/ds1['ro_p']*5.0/6.0*T0/T1/1000.0,label='T_p',color='b',linewidth=lthick)
axs[0,1].plot(ds1['xgrid']/ds1['time'],ds1['pr_n']/ds1['ro_n']*5.0/3.0*T0/T1/1000.0,label='T_n',color='r',linewidth=lthick)
axs[0,1].plot([-1,4],[Tmhd,Tmhd],'k--',label='T (MHD)',linewidth=lthick)
#axs[0,1].set_yscale('log')
axs[0,1].set_xscale('log')
axs[0,1].set_xlim([0.001,4])#axs[0,1].set_xlim([0.01,2.0])
#axs[0,1].legend()
#axs[0,1].set_xlabel('$x_s$')
axs[0,1].set_ylabel('Temperature [kK]')

axs[0,2].plot(ds1['xgrid']/ds1['time'],ds1['ion_loss'],label='P_n',color='k',linewidth=lthick)
#axs[1,1].plot(xs,ds['pr_p'],label='P_p',color='b',linewidth=lthick)
#axs[1,1].plot(xs,ds['pr_p']+ds['pr_n'],label='P_p+P_n',color='g',linewidth=lthick)
axs[0,2].set_xscale('log')
axs[0,2].set_xlim([0.001,4])
#axs[0,0].set_ylim([-0.16,0.16])
#axs[1,1].legend()
#axs[0,2].set_xlabel('$x_s$')
axs[0,2].set_ylabel('loss rate')


axs[0,3].plot(ds1['xgrid']/ds1['time'],ds1['ion'],label='ion',linewidth=lthick)
axs[0,3].plot(ds1['xgrid']/ds1['time'],ds1['rec'],label='rec',linewidth=lthick)
axs[0,3].plot(ds1['xgrid']/ds1['time'],ds1['ion_rad'],label='ion_rad',linewidth=lthick)
axs[0,3].plot(ds1['xgrid']/ds1['time'],ds1['rec_rad'],label='rec_rad',linewidth=lthick)
axs[0,3].set_yscale('log')
axs[0,3].set_xscale('log')
axs[0,3].set_xlim([0.001,4])
#axs[0,3].legend()
#axs[0,3].set_xlabel('$x_s$')
axs[0,3].set_ylabel('$\Gamma$')

##############################################################################################
axs[1,0].plot(ds10['xgrid']/ds10['time'],ds10['vx_n'],label='vx_n',color='r',linewidth=lthick)
axs[1,0].plot(ds10['xgrid']/ds10['time'],ds10['vx_p'],label='vx_p',color='b',linewidth=lthick)
axs[1,0].set_xscale('log')
axs[1,0].set_xlim([0.001,4])
#axs[0,0].set_ylim([-0.16,0.16])
#axs[1,0].legend()
#axs[1,0].set_xlabel('$x/t$')
axs[1,0].set_ylabel('$v_x$')

T1=ds10['pr_p'][-1]/ds10['ro_p'][-1]*5.0/6.0
T1m=dsm['pr_p'][-1]/dsm['ro_p'][-1]*5.0/3.0
Tmhd=dsm['pr_p'][750]/dsm['ro_p'][750]*5.0/3.0*T0/T1m/1000.0
axs[1,1].plot(ds10['xgrid']/ds10['time'],ds10['pr_p']/ds10['ro_p']*5.0/6.0*T0/T1/1000.0,label='T_p',color='b',linewidth=lthick)
axs[1,1].plot(ds10['xgrid']/ds10['time'],ds10['pr_n']/ds10['ro_n']*5.0/3.0*T0/T1/1000.0,label='T_n',color='r',linewidth=lthick)
axs[1,1].plot([-1,4],[Tmhd,Tmhd],'k--',label='T (MHD)',linewidth=lthick)
#axs[0,1].set_yscale('log')
axs[1,1].set_xscale('log')
axs[1,1].set_xlim([0.001,4])#axs[0,1].set_xlim([0.01,2.0])
#axs[1,1].legend()
#axs[1,1].set_xlabel('$x_s$')
axs[1,1].set_ylabel('Temperature [kK]')

axs[1,2].plot(ds10['xgrid']/ds10['time'],ds10['ion_loss'],label='P_n',color='k',linewidth=lthick)
#axs[1,1].plot(xs,ds['pr_p'],label='P_p',color='b',linewidth=lthick)
#axs[1,1].plot(xs,ds['pr_p']+ds['pr_n'],label='P_p+P_n',color='g',linewidth=lthick)
axs[1,2].set_xscale('log')
axs[1,2].set_xlim([0.001,4])
#axs[0,0].set_ylim([-0.16,0.16])
#axs[1,1].legend()
#axs[1,2].set_xlabel('$x_s$')
axs[1,2].set_ylabel('loss rate')


axs[1,3].plot(ds10['xgrid']/ds10['time'],ds10['ion'],label='ion',linewidth=lthick)
axs[1,3].plot(ds10['xgrid']/ds10['time'],ds10['rec'],label='rec',linewidth=lthick)
axs[1,3].plot(ds10['xgrid']/ds10['time'],ds10['ion_rad'],label='ion_rad',linewidth=lthick)
axs[1,3].plot(ds10['xgrid']/ds10['time'],ds10['rec_rad'],label='rec_rad',linewidth=lthick)
axs[1,3].set_yscale('log')
axs[1,3].set_xscale('log')
axs[1,3].set_xlim([0.001,4])
#axs[0,3].legend()
#axs[1,3].set_xlabel('$x_s$')
axs[1,3].set_ylabel('$\Gamma$')

##############################################################################################
axs[2,0].plot(ds100['xgrid']/ds100['time'],ds100['vx_n'],label='vx_n',color='r',linewidth=lthick)
axs[2,0].plot(ds100['xgrid']/ds100['time'],ds100['vx_p'],label='vx_p',color='b',linewidth=lthick)
axs[2,0].set_xscale('log')
axs[2,0].set_xlim([0.001,4])
#axs[0,0].set_ylim([-0.16,0.16])
#axs[1,0].legend()
axs[2,0].set_xlabel('$x/t$')
axs[2,0].set_ylabel('$v_x$')

T1=ds100['pr_p'][-1]/ds100['ro_p'][-1]*5.0/6.0
T1m=dsm['pr_p'][-1]/dsm['ro_p'][-1]*5.0/3.0
Tmhd=dsm['pr_p'][750]/dsm['ro_p'][750]*5.0/3.0*T0/T1m/1000.0
axs[2,1].plot(ds100['xgrid']/ds100['time'],ds100['pr_p']/ds100['ro_p']*5.0/6.0*T0/T1/1000.0,label='T_p',color='b',linewidth=lthick)
axs[2,1].plot(ds100['xgrid']/ds100['time'],ds100['pr_n']/ds100['ro_n']*5.0/3.0*T0/T1/1000.0,label='T_n',color='r',linewidth=lthick)
axs[2,1].plot([-1,4],[Tmhd,Tmhd],'k--',label='T (MHD)',linewidth=lthick)
#axs20,1].set_yscale('log')
axs[2,1].set_xscale('log')
axs[2,1].set_xlim([0.001,4])#axs[0,1].set_xlim([0.01,2.0])
#axs[1,1].legend()
axs[2,1].set_xlabel('$x/t$')
axs[2,1].set_ylabel('Temperature [kK]')

axs[2,2].plot(ds100['xgrid']/ds100['time'],ds100['ion_loss'],label='P_n',color='k',linewidth=lthick)
#axs[1,1].plot(xs,ds['pr_p'],label='P_p',color='b',linewidth=lthick)
#axs[1,1].plot(xs,ds['pr_p']+ds['pr_n'],label='P_p+P_n',color='g',linewidth=lthick)
axs[2,2].set_xscale('log')
axs[2,2].set_xlim([0.001,4])
#axs[0,0].set_ylim([-0.16,0.16])
#axs[1,1].legend()
axs[2,2].set_xlabel('$x/t$')
axs[2,2].set_ylabel('loss rate')


axs[2,3].plot(ds100['xgrid']/ds100['time'],ds100['ion'],label='ion',linewidth=lthick)
axs[2,3].plot(ds100['xgrid']/ds100['time'],ds100['rec'],label='rec',linewidth=lthick)
axs[2,3].plot(ds100['xgrid']/ds100['time'],ds100['ion_rad'],label='ion_rad',linewidth=lthick)
axs[2,3].plot(ds100['xgrid']/ds100['time'],ds100['rec_rad'],label='rec_rad',linewidth=lthick)
axs[2,3].set_yscale('log')
axs[2,3].set_xscale('log')
axs[2,3].set_xlim([0.001,4])
#axs[0,3].legend()
axs[2,3].set_xlabel('$x/t$')
axs[2,3].set_ylabel('$\Gamma$')

plt.savefig('timeseries_plot_upc.png',dpi=300)
