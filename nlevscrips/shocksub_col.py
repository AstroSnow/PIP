#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 09:23:02 2022

@author: snow
"""

import pipreadmods as PIPpy
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'font.size': 16})

#fname='../Data/'
#fname='../simdata/nleveltest_rad/Data/'
#fname='../simdata/nleveltest_rad_thick_st'
fname='../simdata/schemtest_pip/' ;xmin=-0.001;xmax=0.0015
#fname='../simdata/isca_midc_ltest/Data/' ;xmin=0.052;xmax=0.074

ds=PIPpy.pipread(fname,18)
#dsm=PIPpy.pipread('../simdata/MHD_ref/',50)


#############################################################
#Estimate the shock frame
b=np.argmin(np.gradient(ds['vx_p'])) #plasma shock location
vs=ds['xgrid'][b]/ds['time'] #plasma shock velocity
xs=ds['xgrid']/ds['time']-vs #shock grid
bn=np.argmin(np.gradient(ds['vx_n'])) #neutral shock location
vsn=ds['xgrid'][bn]/ds['time'] #neutral shock velocity

#bm=np.argmin(np.gradient(dsm['vx_p'])) #MHD shock location
#vsm=dsm['xgrid'][bm]/dsm['time'][0] #MHD shock velocity
#xsm=dsm['xgrid']/dsm['time'][0]-vsm #MHD shock grid

 

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

fig, axs = plt.subplots(1, 1,dpi=300)
fig.set_size_inches(9.7, 6)
lthick=2.5

#axs[0,0].plot(ds['xgrid']/ds['time'],ds['vx_n']-ds['vx_p'],label='vx_n-vx_p',color='k',linewidth=lthick)
axs.plot(xs,ds['vx_n'],label='vx_n',color='r',linewidth=lthick)
axs.plot(xs,ds['vx_p'],label='vx_p',color='b',linewidth=lthick)
#ax[0,0].set_xscale('log')
axs.fill_between([0.00082,0.002],[-0.3,-0.3], [0.02,0.02], alpha=0.3,facecolor='green')
axs.fill_between([-0.001,-0.0002],[-0.3,-0.3], [0.02,0.02], alpha=0.3,facecolor='green')
axs.text(-0.0008,-0.1,'post-shock',color='k',fontsize=16)
axs.text(0.001,-0.1,'pre-shock',color='k',fontsize=16)
axs.text(0.0,-0.18,'finite-width',color='k',fontsize=16)
#axs.plot([xs[bn],xs[bn]],[-1,1],'r--')
#axs.plot([xs[b],xs[b]],[-1,1],'b--')
axs.set_xlim([xmin,xmax])
axs.set_ylim([-0.20,0.01])
axs.legend()
axs.set_xlabel('$xs$')
axs.set_ylabel('Velocity')

plt.savefig('shocksub_col.png',dpi=300)