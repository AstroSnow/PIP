#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 13:02:45 2022

@author: ben
"""
import pipreadmods as PIPpy
import matplotlib.pyplot as plt
import numpy as np


#fname='../Data/'
#fname='../simdata/nleveltest_rad/Data/'
#fname='../simdata/nleveltest_rad_thick_st'
fname='../simdata/T_6220_long/Data/'; T0=6220.0; etime=30; sname='contextplot.png'
#fname='../simdata/isca_midc/Data/'
#fname='../simdata/isca_midc_ltest/Data/' ; T0=5030.0; etime=40; sname='context_midc.png'
#fname='../simdata/isca_lowc_long/Data/'; T0=5180.0; etime=21; sname='context_lowc.png'


ds=PIPpy.pipread(fname,etime)
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

fig, axs = plt.subplots(2, 2,dpi=300)
fig.set_size_inches(9.7, 6)
lthick=2.5

"""axs[0,1].plot(ds['xgrid']/ds['time'],ds['ion'],label='ion')
axs[0,1].plot(ds['xgrid']/ds['time'],ds['rec'],label='rec')
axs[0,1].plot(ds['xgrid']/ds['time'],ds['ion_rad'],label='ion_rad')
axs[0,1].plot(ds['xgrid']/ds['time'],ds['rec_rad'],label='rec_rad')
axs[0,1].set_yscale('log')
axs[0,1].set_xscale('log')
axs[0,1].set_xlim([0.01,2.0])
axs[0,1].legend()
axs[0,1].set_xlabel('$x/t$')
axs[0,1].set_ylabel('$\Gamma$')"""

axs[0,0].plot(ds['xgrid']/ds['time'],ds['vx_p'],label='vx_p',color='b',linewidth=lthick)
axs[0,0].plot(ds['xgrid']/ds['time'],ds['vx_n'],label='vx_n',color='r',linewidth=lthick)
axs[0,0].plot(dsm['xgrid']/dsm['time'],dsm['vx_p'],label='vx (MHD)',color='k',linewidth=lthick)
axs[0,0].set_xscale('log')
axs[0,0].set_xlim([0.01,2.0])
axs[0,0].text(0.2,-0.05,'rarefaction',color='k')
axs[0,0].text(0.28,-0.07,'wave',color='k')
axs[0,0].arrow(0.4,-0.08,0.3,-0.05,color='g',width=0.002,head_width=0.02,head_length=0.05)
axs[0,0].text(0.18,-0.23,'inflow',color='k')
axs[0,0].arrow(0.25,-0.215,0.0,0.02,color='g',width=0.002,head_width=0.02,head_length=0.01)
axs[0,0].arrow(0.25,-0.235,0.0,-0.02,color='g',width=0.002,head_width=0.02,head_length=0.01)
axs[0,0].text(0.012,-0.1,'slow-mode',color='k')
axs[0,0].text(0.015,-0.12,'shock',color='k')
axs[0,0].arrow(0.032,-0.115,0.03,0.0,color='g',width=0.002,head_width=0.02,head_length=0.01)
#axs[0,0].arrow(0.25,-0.235,0.0,-0.02,color='g',width=0.002,head_width=0.02,head_length=0.01)
axs[0,0].legend()
axs[0,0].set_xlabel('$x/t$')
axs[0,0].set_ylabel('$v_x$')

axs[0,1].plot(ds['xgrid']/ds['time'],ds['by'],label='B_y',color='b',linewidth=lthick)
axs[0,1].plot(dsm['xgrid']/dsm['time'],dsm['by'],label='B_y (MHD)',color='k',linewidth=lthick)
#axs[0,1].set_yscale('log')
axs[0,1].set_xscale('log')
axs[0,1].set_xlim([0.01,2.0])
axs[0,1].legend()
axs[0,1].set_xlabel('$x/t$')
axs[0,1].set_ylabel('$B_y$')

T1=ds['pr_p'][-1]/ds['ro_p'][-1]*5.0/6.0
T1m=dsm['pr_p'][-1]/dsm['ro_p'][-1]*5.0/3.0

axs[1,0].plot(ds['xgrid']/ds['time'],ds['pr_p']/ds['ro_p']*5.0/6.0*T0/T1,label='T_p',color='b',linewidth=lthick)
axs[1,0].plot(ds['xgrid']/ds['time'],ds['pr_n']/ds['ro_n']*5.0/3.0*T0/T1,label='T_n',color='r',linewidth=lthick)
axs[1,0].plot(dsm['xgrid']/dsm['time'],dsm['pr_p']/dsm['ro_p']*5.0/3.0*T0/T1m,label='T_p (MHD)',color='k',linewidth=lthick)
axs[1,0].legend()
axs[1,0].set_xlim([0.01,2.0])
axs[1,0].set_xscale('log')
axs[1,0].set_xlabel('$x/t$')
axs[1,0].set_ylabel('Temperature [K]')

"""axs[1,1].plot(ds['xgrid']/ds['time'],ds['nexcite1'],label='n1')
axs[1,1].plot(ds['xgrid']/ds['time'],ds['nexcite2'],label='n2')
axs[1,1].plot(ds['xgrid']/ds['time'],ds['nexcite3'],label='n3')
axs[1,1].plot(ds['xgrid']/ds['time'],ds['nexcite4'],label='n4')
axs[1,1].plot(ds['xgrid']/ds['time'],ds['nexcite5'],label='n5')
axs[1,1].plot(ds['xgrid']/ds['time'],ds['nexcite6'],label='c')
axs[1,1].legend()
axs[1,1].set_xscale('log')
axs[1,1].set_yscale('log')
axs[1,1].set_xlim([0.01,2.0])
axs[1,1].set_xlabel('$x/t$')
axs[1,1].set_ylabel('Density')"""

axs[1,1].plot(ds['xgrid']/ds['time'],ds['ro_p'],label='ro_p',color='b',linewidth=lthick)
axs[1,1].plot(ds['xgrid']/ds['time'],ds['ro_n'],label='ro_n',color='r',linewidth=lthick)
axs[1,1].plot(ds['xgrid']/ds['time'],ds['ro_n']+ds['ro_p'],label='ro_n+ro_p',color='g',linewidth=lthick)
axs[1,1].plot(dsm['xgrid']/dsm['time'],dsm['ro_p'],label='ro_p (MHD)',color='k',linewidth=lthick)
axs[1,1].legend()
axs[1,1].set_xscale('log')
#axs[1,1].set_yscale('log')
axs[1,1].set_xlim([0.01,2.0])
axs[1,1].set_xlabel('$x/t$')
axs[1,1].set_ylabel('Density')

#plt.savefig('contextplot.png',dpi=300)
plt.savefig(sname,dpi=300)

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
