#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 24 08:00:40 2022

@author: ben
"""
import pipreadmods as PIPpy
import matplotlib.pyplot as plt
import numpy as np

#fname='../simdata/col_t10000/' ;sname='col_10000' #xi_n=0.94787198761607150
fname='../simdata/rad_thin_t10000/' ;sname='rad_thin_10000' #xi_n=0.94787198863788269
#fname='../simdata/rad_thick_t10000/' ;sname='rad_thick_10000' #xi_n=0.94787198863788269

ds=PIPpy.pipread(fname,10)

dsm=PIPpy.pipread('../simdata/MHD_t10000/',10) # Reference MHD simulation

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

fig, axs = plt.subplots(2, 2)
fig.set_size_inches(12.4, 7)
fig.tight_layout(pad=2.8)
lthick=2.5

axs[0,0].plot(ds['xgrid']/ds['time'],ds['ion'],label='ion',linewidth=lthick,color='orange')
axs[0,0].plot(ds['xgrid']/ds['time'],ds['rec'],label='rec',linewidth=lthick,color='green')
if 'ion_rad' in ds:
    axs[0,0].plot(ds['xgrid']/ds['time'],ds['ion_rad'],label='ion_rad',linewidth=lthick,color='c')
    axs[0,0].plot(ds['xgrid']/ds['time'],ds['rec_rad'],label='rec_rad',linewidth=lthick,color='y')
axs[0,0].set_yscale('log')
axs[0,0].set_xscale('log')
axs[0,0].legend()
axs[0,0].set_xlim([0.001,1.5])
axs[0,0].set_xlabel('x/t')
axs[0,0].set_ylabel('$\Gamma$')

axs[0,1].plot(ds['xgrid']/ds['time'],ds['vx_p'],label='$v_{xp}$',linewidth=lthick,color='b')
axs[0,1].plot(ds['xgrid']/ds['time'],ds['vx_n'],label='$v_{xn}$',linewidth=lthick,color='r')
axs[0,1].plot(dsm['xgrid']/dsm['time'],dsm['vx_p'],label='$v_{x}$ (MHD)',linewidth=lthick,color='k')
axs[0,1].set_xscale('log')
axs[0,1].legend()
axs[0,1].set_xlim([0.001,1.5])
axs[0,1].set_xlabel('x/t')
axs[0,1].set_ylabel('$v_{x}$')

axs[1,0].plot(ds['xgrid']/ds['time'],ds['pr_p']/ds['ro_p']*5.0/6.0,label='$T_p$',linewidth=lthick,color='b')
axs[1,0].plot(ds['xgrid']/ds['time'],ds['pr_n']/ds['ro_n']*5.0/3.0,label='$T_n$',linewidth=lthick,color='r')
axs[1,0].plot(dsm['xgrid']/dsm['time'],dsm['pr_p']/dsm['ro_p']*5.0/3.0,label='$T$ (MHD)',linewidth=lthick,color='k')
axs[1,0].legend()
axs[1,0].set_xscale('log')
axs[1,0].set_xlim([0.001,1.5])
axs[1,0].set_xlabel('x/t')
axs[1,0].set_ylabel('Temperature')

axs[1,1].plot(ds['xgrid']/ds['time'],ds['nexcite1'],label='n1',linewidth=lthick,color='r')
axs[1,1].plot(ds['xgrid']/ds['time'],ds['nexcite2'],label='n2',linewidth=lthick)
axs[1,1].plot(ds['xgrid']/ds['time'],ds['nexcite3'],label='n3',linewidth=lthick)
axs[1,1].plot(ds['xgrid']/ds['time'],ds['nexcite4'],label='n4',linewidth=lthick)
axs[1,1].plot(ds['xgrid']/ds['time'],ds['nexcite5'],label='n5',linewidth=lthick)
axs[1,1].plot(ds['xgrid']/ds['time'],ds['nexcite6'],label='p',linewidth=lthick,color='b')
axs[1,1].legend()
axs[1,1].set_xscale('log')
axs[1,1].set_yscale('log')
axs[1,1].set_xlim([0.001,1.5])
axs[1,1].set_xlabel('x/t')
axs[1,1].set_ylabel('Level populations')

plt.savefig(''.join(['simlevels_',sname]),dpi=300)