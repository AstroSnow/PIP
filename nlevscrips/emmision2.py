#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 24 09:58:33 2022

@author: ben
"""

import pipreadmods as PIPpy
import matplotlib.pyplot as plt
import numpy as np

#fname='../simdata/col_t10000/' ;sname='col_10000' #xi_n=0.94787198761607150
#fname='../simdata/rad_thin_t10000/' ;sname='rad_thin_10000' #xi_n=0.94787198863788269
#fname='../simdata/rad_thick_t10000/' ;sname='rad_thick_10000' #xi_n=0.94787198863788269
fname='../simdata/T_6220_long/Data/'

ds=PIPpy.pipread(fname,30,exrates=1)

fig, axs = plt.subplots(2, 1)
fig.set_size_inches(12.4, 7)
fig.tight_layout(pad=2.8)
lthick=2.5

"""axs[0].plot(ds['xgrid']/ds['time'],ds['colrat'][:,2,1]*ds['nexcite2'],label='$\Gamma_{2,1} n_2$',linewidth=lthick)
axs[0].plot(ds['xgrid']/ds['time'],ds['colrat'][:,3,1]*ds['nexcite3'],label='$\Gamma_{3,1} n_3$',linewidth=lthick)
axs[0].plot(ds['xgrid']/ds['time'],ds['colrat'][:,4,1]*ds['nexcite4'],label='$\Gamma_{4,1} n_4$',linewidth=lthick)
axs[0].plot(ds['xgrid']/ds['time'],ds['colrat'][:,5,1]*ds['nexcite5'],label='$\Gamma_{5,1} n_5$',linewidth=lthick)
axs[0].set_yscale('log')
axs[0].set_xscale('log')
axs[0].legend()
axs[0].set_xlim([0.001,1.5])
axs[0].set_xlabel('x/t')
axs[0].set_ylabel('Lyman')

axs[1].plot(ds['xgrid']/ds['time'],ds['colrat'][:,3,2]*ds['nexcite3'],label='$\Gamma_{3,2} n_3$',linewidth=lthick)
axs[1].plot(ds['xgrid']/ds['time'],ds['colrat'][:,4,2]*ds['nexcite4'],label='$\Gamma_{4,2} n_4$',linewidth=lthick)
axs[1].plot(ds['xgrid']/ds['time'],ds['colrat'][:,5,2]*ds['nexcite5'],label='$\Gamma_{5,2} n_5$',linewidth=lthick)
axs[1].set_xscale('log')
axs[1].set_yscale('log')
axs[1].legend()
axs[1].set_xlim([0.001,1.5])
axs[1].set_xlabel('x/t')
axs[1].set_ylabel('Balmer')"""

####################################

axs[0].plot(ds['xgrid']/ds['time'],(ds['colrat'][:,2,1]+ds['radrat'][:,2,1])*ds['nexcite2'],label='$\Gamma_{2,1} n_2$',linewidth=lthick)
axs[0].plot(ds['xgrid']/ds['time'],(ds['colrat'][:,3,1]+ds['radrat'][:,3,1])*ds['nexcite3'],label='$\Gamma_{3,1} n_3$',linewidth=lthick)
axs[0].plot(ds['xgrid']/ds['time'],(ds['colrat'][:,4,1]+ds['radrat'][:,4,1])*ds['nexcite4'],label='$\Gamma_{4,1} n_4$',linewidth=lthick)
axs[0].plot(ds['xgrid']/ds['time'],(ds['colrat'][:,5,1]+ds['radrat'][:,5,1])*ds['nexcite5'],label='$\Gamma_{5,1} n_5$',linewidth=lthick)
axs[0].set_yscale('log')
axs[0].set_xscale('log')
axs[0].legend()
axs[0].set_xlim([0.075,0.08])
axs[0].set_xlabel('x/t')
axs[0].set_ylabel('Lyman')

axs[1].plot(ds['xgrid']/ds['time'],(ds['colrat'][:,3,2]+ds['radrat'][:,3,2])*ds['nexcite3'],label='$\Gamma_{3,2} n_3$',linewidth=lthick)
axs[1].plot(ds['xgrid']/ds['time'],(ds['colrat'][:,4,2]+ds['radrat'][:,4,2])*ds['nexcite4'],label='$\Gamma_{4,2} n_4$',linewidth=lthick)
axs[1].plot(ds['xgrid']/ds['time'],(ds['colrat'][:,5,2]+ds['radrat'][:,5,2])*ds['nexcite5'],label='$\Gamma_{5,2} n_5$',linewidth=lthick)
axs[1].set_xscale('log')
axs[1].set_yscale('log')
axs[1].legend()
axs[1].set_xlim([0.075,0.08])
axs[1].set_xlabel('x/t')
axs[1].set_ylabel('Balmer')


#plt.savefig('emmision2_t10000_uc.png', dpi=300)