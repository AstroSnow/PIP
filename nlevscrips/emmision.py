#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 24 09:58:33 2022

@author: ben
"""

import pipreadmods as PIPpy
import matplotlib.pyplot as plt
import numpy as np

fname='../simdata/col_t10000/' ;sname='col_10000' #xi_n=0.94787198761607150
#fname='../simdata/rad_thin_t10000/' ;sname='rad_thin_10000' #xi_n=0.94787198863788269
#fname='../simdata/rad_thick_t10000/' ;sname='rad_thick_10000' #xi_n=0.94787198863788269

ds=PIPpy.pipread(fname,10,exrates=1)

fig, axs = plt.subplots(3, 2)
fig.set_size_inches(12.4, 7)
fig.tight_layout(pad=2.8)
lthick=2.5

axs[0,0].plot(ds['xgrid']/ds['time'],ds['colrat'][:,2,1],label='col21',linewidth=lthick)
axs[0,0].plot(ds['xgrid']/ds['time'],ds['colrat'][:,3,1],label='col31',linewidth=lthick)
axs[0,0].plot(ds['xgrid']/ds['time'],ds['colrat'][:,4,1],label='col41',linewidth=lthick)
axs[0,0].plot(ds['xgrid']/ds['time'],ds['colrat'][:,5,1],label='col51',linewidth=lthick)
axs[0,0].set_yscale('log')
axs[0,0].set_xscale('log')
axs[0,0].legend()
axs[0,0].set_xlim([0.001,1.5])
axs[0,0].set_xlabel('x/t')
axs[0,0].set_ylabel('Lyman')

axs[0,1].plot(ds['xgrid']/ds['time'],ds['colrat'][:,2,1],label='col32',linewidth=lthick)
axs[0,1].plot(ds['xgrid']/ds['time'],ds['colrat'][:,3,1],label='col42',linewidth=lthick)
axs[0,1].plot(ds['xgrid']/ds['time'],ds['colrat'][:,4,1],label='col52',linewidth=lthick)
axs[0,1].set_xscale('log')
axs[0,1].set_yscale('log')
axs[0,1].legend()
axs[0,1].set_xlim([0.001,1.5])
axs[0,1].set_xlabel('x/t')
axs[0,1].set_ylabel('Balmer')

###############################################################################
fname='../simdata/rad_thin_t10000/'
ds=PIPpy.pipread(fname,10,exrates=1)

axs[1,0].plot(ds['xgrid']/ds['time'],ds['colrat'][:,2,1]+ds['radrat'][:,2,1],label='col21+rad21',linewidth=lthick)
axs[1,0].plot(ds['xgrid']/ds['time'],ds['colrat'][:,3,1]+ds['radrat'][:,3,1],label='col31+rad31',linewidth=lthick)
axs[1,0].plot(ds['xgrid']/ds['time'],ds['colrat'][:,4,1]+ds['radrat'][:,4,1],label='col41+rad41',linewidth=lthick)
axs[1,0].plot(ds['xgrid']/ds['time'],ds['colrat'][:,5,1]+ds['radrat'][:,5,1],label='col51+rad51',linewidth=lthick)
axs[1,0].set_yscale('log')
axs[1,0].set_xscale('log')
axs[1,0].legend()
axs[1,0].set_xlim([0.001,1.5])
axs[1,0].set_xlabel('x/t')
axs[1,0].set_ylabel('Lyman')

axs[1,1].plot(ds['xgrid']/ds['time'],ds['colrat'][:,2,1]+ds['radrat'][:,3,2],label='col32+rad32',linewidth=lthick)
axs[1,1].plot(ds['xgrid']/ds['time'],ds['colrat'][:,3,1]+ds['radrat'][:,4,2],label='col42+rad42',linewidth=lthick)
axs[1,1].plot(ds['xgrid']/ds['time'],ds['colrat'][:,4,1]+ds['radrat'][:,5,2],label='col52+rad52',linewidth=lthick)
axs[1,1].set_xscale('log')
axs[1,1].set_yscale('log')
axs[1,1].legend()
axs[1,1].set_xlim([0.001,1.5])
axs[1,1].set_xlabel('x/t')
axs[1,1].set_ylabel('Balmer')

###############################################################################
fname='../simdata/rad_thick_t10000/'
ds=PIPpy.pipread(fname,10,exrates=1)

axs[2,0].plot(ds['xgrid']/ds['time'],ds['colrat'][:,2,1]+ds['radrat'][:,2,1],label='col21+rad21',linewidth=lthick)
axs[2,0].plot(ds['xgrid']/ds['time'],ds['colrat'][:,3,1]+ds['radrat'][:,3,1],label='col31+rad31',linewidth=lthick)
axs[2,0].plot(ds['xgrid']/ds['time'],ds['colrat'][:,4,1]+ds['radrat'][:,4,1],label='col41+rad41',linewidth=lthick)
axs[2,0].plot(ds['xgrid']/ds['time'],ds['colrat'][:,5,1]+ds['radrat'][:,5,1],label='col51+rad51',linewidth=lthick)
axs[2,0].set_yscale('log')
axs[2,0].set_xscale('log')
axs[2,0].legend()
axs[2,0].set_xlim([0.001,1.5])
axs[2,0].set_xlabel('x/t')
axs[2,0].set_ylabel('Lyman')

axs[2,1].plot(ds['xgrid']/ds['time'],ds['colrat'][:,2,1]+ds['radrat'][:,3,2],label='col32+rad32',linewidth=lthick)
axs[2,1].plot(ds['xgrid']/ds['time'],ds['colrat'][:,3,1]+ds['radrat'][:,4,2],label='col42+rad42',linewidth=lthick)
axs[2,1].plot(ds['xgrid']/ds['time'],ds['colrat'][:,4,1]+ds['radrat'][:,5,2],label='col52+rad52',linewidth=lthick)
axs[2,1].set_xscale('log')
axs[2,1].set_yscale('log')
axs[2,1].legend()
axs[2,1].set_xlim([0.001,1.5])
axs[2,1].set_xlabel('x/t')
axs[2,1].set_ylabel('Balmer')

plt.savefig('emmision_t10000.png', dpi=300)