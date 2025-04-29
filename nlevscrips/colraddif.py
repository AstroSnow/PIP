#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 11:23:55 2022

@author: ben
"""
import pipreadmods as PIPpy
import matplotlib.pyplot as plt
import numpy as np

fname='../simdata/col_t10000/' ;sname='col_10000' #xi_n=0.94787198761607150
fname2='../simdata/rad_thin_t10000/' ;sname='rad_thin_10000' #xi_n=0.94787198863788269
#fname='../simdata/rad_thick_t10000/' ;sname='rad_thick_10000' #xi_n=0.94787198863788269

ds=PIPpy.pipread(fname,10)
ds2=PIPpy.pipread(fname2,10)

fig, axs = plt.subplots(1,1)
fig.set_size_inches(9.7, 6)

lthick=2.5

var='ion_loss'

axs.plot(ds['xgrid']/ds['time'],(ds[var]-ds2[var])/ds[var],linewidth=lthick)
#axs.set_yscale('log')
#axs.plot(xs,abs(ds['vx_p']-vs)/va,'b',label='Alfven',linewidth=lthick)
#axs.plot(xs,abs(ds['vx_n']-vsn)/csn,'r',label='Sonic',linewidth=lthick)

axs.set_xlabel('$x/t$')
#axs.set_ylabel('$\\rho_p \Gamma_{rec} - \\rho_n \Gamma_{ion}$')
axs.set_xlim([0.02,0.04])
#axs.legend()

#plt.savefig(''.join(['masschange_',sname]),dpi=300)