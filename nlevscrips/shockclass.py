#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 12:04:27 2022

@author: ben
"""
import pipreadmods as PIPpy
import matplotlib.pyplot as plt
import numpy as np


#fname='../Data/'
#fname='../simdata/nleveltest_rad/Data/'
#fname='../simdata/nleveltest_rad_thick_st'
fname='../simdata/T_6220_long/Data/' ;xmin=0.075;xmax=0.08 ; T0=6220.0
#fname='../simdata/isca_midc_ltest/Data/' ;xmin=0.052;xmax=0.074

ds=PIPpy.pipread(fname,30)
dsm=PIPpy.pipread('../simdata/MHD_ref/',50)

nntot=ds['nexcite1']+ds['nexcite2']+ds['nexcite3']+ds['nexcite4']+ds['nexcite5']

ndif=ds['ro_n']-nntot


#############################################################
#Estimate the shock frame
b=np.argmin(np.gradient(ds['vx_p'])) #plasma shock location
vs=ds['xgrid'][b]/ds['time'] #plasma shock velocity
xs=ds['xgrid']/ds['time']-vs #shock grid
bn=np.argmin(np.gradient(ds['vx_n'])) #neutral shock location
vsn=ds['xgrid'][bn]/ds['time'] #neutral shock velocity

bm=np.argmin(np.gradient(dsm['vx_p'])) #MHD shock location
vsm=dsm['xgrid'][bm]/dsm['time'] #MHD shock velocity
xsm=dsm['xgrid']/dsm['time']-vsm #MHD shock grid

############################################################
#Wave speeds
va=ds['bx']/np.sqrt(ds['ro_p']+ds['ro_n'])
vap=ds['bx']/np.sqrt(ds['ro_p'])

csn=np.sqrt(5.0/3.0*ds['pr_n']/ds['ro_n'])
csp=np.sqrt(5.0/3.0*ds['pr_p']/ds['ro_p'])

fig, axs = plt.subplots(2, 2,dpi=300)
fig.set_size_inches(9.7, 6)
lthick=2.5

#axs[0,0].plot(ds['xgrid']/ds['time'],ds['vx_n']-ds['vx_p'],label='vx_n-vx_p',color='k',linewidth=lthick)
axs[0,0].plot(xs,(ds['vx_n']-vsn)/csn,label='Sonic',color='r',linewidth=lthick)
axs[0,0].plot(xs,(ds['vx_p']-vs)/va,label='Bulk Alfven',color='g',linewidth=lthick)
axs[0,0].plot(xs,(ds['vx_p']-vs)/vap,label='Plasma Alfven',color='b',linewidth=lthick)
#axs[0,0].set_xscale('log')
axs[0,0].set_xlim([xmin-vs,xmax-vs])
#axs[0,0].set_ylim([-0.16,0.16])
axs[0,0].legend()
axs[0,0].set_xlabel('$x_s$')
axs[0,0].set_ylabel('Mach number')

#axs[0,1].plot(xs,ds['ro_n']/(ds['ro_n']+ds['ro_p']),color='k',linewidth=lthick)
axs[0,1].plot(xs,ds['ro_p'],color='b',linewidth=lthick)
axs[0,1].plot(xs,ds['ro_n'],color='r',linewidth=lthick)
#axs[0,1].plot(xs,(ds['vx_p']-vs)/va,label='vx_p-vs',color='b',linewidth=lthick)
#axs[0,1].set_yscale('log')
axs[0,1].set_xlim([xmin-vs,xmax-vs])
#axs[0,1].set_ylim([0.1,1.0])
#axs[0,1].legend()
axs[0,1].set_xlabel('$x_s$')
axs[0,1].set_ylabel('$density$')

axs[1,0].plot(xs,csn,label='cs_n',color='r',linewidth=lthick)
axs[1,0].plot(xs,csp,label='cs_p',color='y',linewidth=lthick)
axs[1,0].plot(xs,va,label='vA',color='g',linewidth=lthick)
axs[1,0].plot(xs,vap,label='vA_p',color='b',linewidth=lthick)
#axs[0,0].set_xscale('log')
axs[1,0].set_xlim([xmin-vs,xmax-vs])
#axs[0,0].set_ylim([-0.16,0.16])
axs[1,0].legend()
axs[1,0].set_xlabel('$x_s$')
axs[1,0].set_ylabel('Wave speeds')

"""axs[1,1].plot(xs,ds['pr_n'],label='P_n',color='r',linewidth=lthick)
axs[1,1].plot(xs,ds['pr_p'],label='P_p',color='b',linewidth=lthick)
axs[1,1].plot(xs,ds['pr_p']+ds['pr_n'],label='P_p+P_n',color='g',linewidth=lthick)
#axs[0,0].set_xscale('log')
axs[1,1].set_xlim([xmin-vs,xmax-vs])
#axs[0,0].set_ylim([-0.16,0.16])
axs[1,1].legend()
axs[1,1].set_xlabel('$x_s$')
axs[1,1].set_ylabel('Pressure')"""

axs[1,1].plot(xs,ds['ion_loss'],label='P_n',color='k',linewidth=lthick)
#axs[1,1].plot(xs,ds['pr_p'],label='P_p',color='b',linewidth=lthick)
#axs[1,1].plot(xs,ds['pr_p']+ds['pr_n'],label='P_p+P_n',color='g',linewidth=lthick)
#axs[0,0].set_xscale('log')
axs[1,1].set_xlim([xmin-vs,xmax-vs])
#axs[0,0].set_ylim([-0.16,0.16])
#axs[1,1].legend()
axs[1,1].set_xlabel('$x_s$')
axs[1,1].set_ylabel('loss rate')

plt.savefig('shockclass_plot_upc.png',dpi=300)
