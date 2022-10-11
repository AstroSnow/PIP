#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 24 08:00:40 2022

@author: ben
"""
import pipreadmods as PIPpy
import matplotlib.pyplot as plt
import numpy as np

#fname='../simdata/col_t10' ;sname='col_et' #xi_n=0.94787198761607150
#fname='../simdata/rad_thin_t10' ;sname='rad_thin_et' #xi_n=0.94787198863788269
fname='../simdata/rad_thick_t10' ;sname='rad_thick_et' #xi_n=0.94787198863788269

ds=PIPpy.pipread(''.join([fname,'/']),1)

fig, axs = plt.subplots(4, 3)
fig.set_size_inches(12.4, 14)
fig.tight_layout(pad=2.8)
lthick=2.5

axs[0,0].plot(ds['xgrid']/ds['time'],ds['vx_p'],label='$v_{xp}$',linewidth=lthick)
axs[0,0].plot(ds['xgrid']/ds['time'],ds['vx_n'],label='$v_{xn}$',linewidth=lthick)
axs[0,0].set_xscale('log')
axs[0,0].legend()
axs[0,0].set_xlim([0.001,15])
axs[0,0].set_xlabel('x/t')
axs[0,0].set_ylabel('$v_{x}$')

axs[1,0].plot(ds['xgrid']/ds['time'],ds['pr_p']/ds['ro_p']*5.0/6.0,label='$T_p$',linewidth=lthick)
axs[1,0].plot(ds['xgrid']/ds['time'],ds['pr_n']/ds['ro_n']*5.0/3.0,label='$T_n$',linewidth=lthick)
axs[1,0].legend()
axs[1,0].set_xscale('log')
axs[1,0].set_xlim([0.001,15])
axs[1,0].set_xlabel('x/t')
axs[1,0].set_ylabel('Temperature')

axs[2,0].plot(ds['xgrid']/ds['time'],ds['ion'],label='ion',linewidth=lthick)
axs[2,0].plot(ds['xgrid']/ds['time'],ds['rec'],label='rec',linewidth=lthick)
if 'ion_rad' in ds:
    axs[2,0].plot(ds['xgrid']/ds['time'],ds['ion_rad'],label='ion_rad',linewidth=lthick)
    axs[2,0].plot(ds['xgrid']/ds['time'],ds['rec_rad'],label='rec_rad',linewidth=lthick)
axs[2,0].set_yscale('log')
axs[2,0].set_xscale('log')
axs[2,0].legend()
axs[2,0].set_xlim([0.001,15])
axs[2,0].set_xlabel('x/t')
axs[2,0].set_ylabel('$\Gamma$')

axs[3,0].plot(ds['xgrid']/ds['time'],ds['nexcite1'],label='n1',linewidth=lthick)
axs[3,0].plot(ds['xgrid']/ds['time'],ds['nexcite2'],label='n2',linewidth=lthick)
axs[3,0].plot(ds['xgrid']/ds['time'],ds['nexcite3'],label='n3',linewidth=lthick)
axs[3,0].plot(ds['xgrid']/ds['time'],ds['nexcite4'],label='n4',linewidth=lthick)
axs[3,0].plot(ds['xgrid']/ds['time'],ds['nexcite5'],label='n5',linewidth=lthick)
axs[3,0].plot(ds['xgrid']/ds['time'],ds['nexcite6'],label='p',linewidth=lthick)
axs[3,0].legend()
axs[3,0].set_xscale('log')
axs[3,0].set_yscale('log')
axs[3,0].set_xlim([0.001,15])
axs[3,0].set_xlabel('x/t')
axs[3,0].set_ylabel('Level populations')

###############################################################################

ds=PIPpy.pipread(''.join([fname,'0/']),1)

axs[0,1].plot(ds['xgrid']/ds['time'],ds['vx_p'],label='$v_{xp}$',linewidth=lthick)
axs[0,1].plot(ds['xgrid']/ds['time'],ds['vx_n'],label='$v_{xn}$',linewidth=lthick)
axs[0,1].set_xscale('log')
axs[0,1].legend()
axs[0,1].set_xlim([0.001,15])
axs[0,1].set_xlabel('x/t')
axs[0,1].set_ylabel('$v_{x}$')

axs[1,1].plot(ds['xgrid']/ds['time'],ds['pr_p']/ds['ro_p']*5.0/6.0,label='$T_p$',linewidth=lthick)
axs[1,1].plot(ds['xgrid']/ds['time'],ds['pr_n']/ds['ro_n']*5.0/3.0,label='$T_n$',linewidth=lthick)
axs[1,1].legend()
axs[1,1].set_xscale('log')
axs[1,1].set_xlim([0.001,15])
axs[1,1].set_xlabel('x/t')
axs[1,1].set_ylabel('Temperature')

axs[2,1].plot(ds['xgrid']/ds['time'],ds['ion'],label='ion',linewidth=lthick)
axs[2,1].plot(ds['xgrid']/ds['time'],ds['rec'],label='rec',linewidth=lthick)
if 'ion_rad' in ds:
    axs[2,1].plot(ds['xgrid']/ds['time'],ds['ion_rad'],label='ion_rad',linewidth=lthick)
    axs[2,1].plot(ds['xgrid']/ds['time'],ds['rec_rad'],label='rec_rad',linewidth=lthick)
axs[2,1].set_yscale('log')
axs[2,1].set_xscale('log')
axs[2,1].legend()
axs[2,1].set_xlim([0.001,15])
axs[2,1].set_xlabel('x/t')
axs[2,1].set_ylabel('$\Gamma$')

axs[3,1].plot(ds['xgrid']/ds['time'],ds['nexcite1'],label='n1',linewidth=lthick)
axs[3,1].plot(ds['xgrid']/ds['time'],ds['nexcite2'],label='n2',linewidth=lthick)
axs[3,1].plot(ds['xgrid']/ds['time'],ds['nexcite3'],label='n3',linewidth=lthick)
axs[3,1].plot(ds['xgrid']/ds['time'],ds['nexcite4'],label='n4',linewidth=lthick)
axs[3,1].plot(ds['xgrid']/ds['time'],ds['nexcite5'],label='n5',linewidth=lthick)
axs[3,1].plot(ds['xgrid']/ds['time'],ds['nexcite6'],label='p',linewidth=lthick)
axs[3,1].legend()
axs[3,1].set_xscale('log')
axs[3,1].set_yscale('log')
axs[3,1].set_xlim([0.001,15])
axs[3,1].set_xlabel('x/t')
axs[3,1].set_ylabel('Level populations')

###############################################################################

ds=PIPpy.pipread(''.join([fname,'0/']),10)

axs[0,2].plot(ds['xgrid']/ds['time'],ds['vx_p'],label='$v_{xp}$',linewidth=lthick)
axs[0,2].plot(ds['xgrid']/ds['time'],ds['vx_n'],label='$v_{xn}$',linewidth=lthick)
axs[0,2].set_xscale('log')
axs[0,2].legend()
axs[0,2].set_xlim([0.001,15])
axs[0,2].set_xlabel('x/t')
axs[0,2].set_ylabel('$v_{x}$')

axs[1,2].plot(ds['xgrid']/ds['time'],ds['pr_p']/ds['ro_p']*5.0/6.0,label='$T_p$',linewidth=lthick)
axs[1,2].plot(ds['xgrid']/ds['time'],ds['pr_n']/ds['ro_n']*5.0/3.0,label='$T_n$',linewidth=lthick)
axs[1,2].legend()
axs[1,2].set_xscale('log')
axs[1,2].set_xlim([0.001,15])
axs[1,2].set_xlabel('x/t')
axs[1,2].set_ylabel('Temperature')

axs[2,2].plot(ds['xgrid']/ds['time'],ds['ion'],label='ion',linewidth=lthick)
axs[2,2].plot(ds['xgrid']/ds['time'],ds['rec'],label='rec',linewidth=lthick)
if 'ion_rad' in ds:
    axs[2,2].plot(ds['xgrid']/ds['time'],ds['ion_rad'],label='ion_rad',linewidth=lthick)
    axs[2,2].plot(ds['xgrid']/ds['time'],ds['rec_rad'],label='rec_rad',linewidth=lthick)
axs[2,2].set_yscale('log')
axs[2,2].set_xscale('log')
axs[2,2].legend()
axs[2,2].set_xlim([0.001,15])
axs[2,2].set_xlabel('x/t')
axs[2,2].set_ylabel('$\Gamma$')

axs[3,2].plot(ds['xgrid']/ds['time'],ds['nexcite1'],label='n1',linewidth=lthick)
axs[3,2].plot(ds['xgrid']/ds['time'],ds['nexcite2'],label='n2',linewidth=lthick)
axs[3,2].plot(ds['xgrid']/ds['time'],ds['nexcite3'],label='n3',linewidth=lthick)
axs[3,2].plot(ds['xgrid']/ds['time'],ds['nexcite4'],label='n4',linewidth=lthick)
axs[3,2].plot(ds['xgrid']/ds['time'],ds['nexcite5'],label='n5',linewidth=lthick)
axs[3,2].plot(ds['xgrid']/ds['time'],ds['nexcite6'],label='p',linewidth=lthick)
axs[3,2].legend()
axs[3,2].set_xscale('log')
axs[3,2].set_yscale('log')
axs[3,2].set_xlim([0.001,15])
axs[3,2].set_xlabel('x/t')
axs[3,2].set_ylabel('Level populations')

plt.savefig(''.join(['timeevo_',sname]),dpi=300)