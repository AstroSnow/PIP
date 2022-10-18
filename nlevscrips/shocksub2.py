#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 13:58:39 2022

@author: ben
"""
import pipreadmods as PIPpy
import matplotlib.pyplot as plt
import numpy as np


#fname='../Data/'
#fname='../simdata/nleveltest_rad/Data/'
#fname='../simdata/nleveltest_rad_thick_st'
#fname='../simdata/T_6220_long/Data/' ;xmin=0.075;xmax=0.08 ;sname='shocksub2_plot_upc.png'; etime=30;T0=6220
#fname='../simdata/isca_midc_ltest/Data/' ;xmin=0.052;xmax=0.074 ;sname='shocksub2_plot_midc.png'; etime=40; T0=5030
fname='../simdata/isca_lowc_long/Data/' ;xmin=0.04;xmax=0.1 ;sname='shocksub2_plot_lowc.png'; etime=21; T0=5180

ds=PIPpy.pipread(fname,etime)
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

#axs[0,0].plot(ds['xgrid']/ds['time'],ds['vx_n']-ds['vx_p'],label='vx_n-vx_p',color='k',linewidth=lthick)
axs[0,0].plot(xs,ds['vx_p'],label='vx_p',color='b',linewidth=lthick)
axs[0,0].plot(xs,ds['vx_n'],label='vx_n',color='r',linewidth=lthick)
#axs[0,0].set_xscale('log')
axs[0,0].set_xlim([xmin-vs,xmax-vs])
#axs[0,0].set_ylim([-0.16,0.16])
axs[0,0].legend()
axs[0,0].set_xlabel('$x_s$')
axs[0,0].set_ylabel('$v_\perp$')

axs[0,1].plot(xs,(ds['vx_n']-vsn)/csn,label='Sonic',color='r',linewidth=lthick)
axs[0,1].plot(xs,(ds['vx_p']-vs)/va,label='Bulk Alfven',color='g',linewidth=lthick)
axs[0,1].plot(xs,(ds['vx_p']-vs)/vap,label='Plasma Alfven',color='b',linewidth=lthick)
#axs[0,0].set_xscale('log')
axs[0,1].set_xlim([xmin-vs,xmax-vs])
#axs[0,0].set_ylim([-0.16,0.16])
axs[0,1].legend()
axs[0,1].set_xlabel('$x_s$')
axs[0,1].set_ylabel('Mach number')

axs[1,0].plot(xs,ds['ion'],label='ion',linewidth=lthick)
axs[1,0].plot(xs,ds['rec'],label='rec',linewidth=lthick)
axs[1,0].plot(xs,ds['ion_rad'],label='ion_rad',linewidth=lthick)
axs[1,0].plot(xs,ds['rec_rad'],label='rec_rad',linewidth=lthick)
axs[1,0].set_yscale('log')
#axs[1,0].set_xscale('log')
axs[1,0].set_xlim([xmin-vs,xmax-vs])
axs[1,0].legend()
axs[1,0].set_xlabel('$x_s$')
axs[1,0].set_ylabel('$\Gamma$')

axs[1,1].plot(xs,ds['nexcite1'],label='n1',linewidth=lthick)
axs[1,1].plot(xs,ds['nexcite2'],label='n2',linewidth=lthick)
axs[1,1].plot(xs,ds['nexcite3'],label='n3',linewidth=lthick)
axs[1,1].plot(xs,ds['nexcite4'],label='n4',linewidth=lthick)
axs[1,1].plot(xs,ds['nexcite5'],label='n5',linewidth=lthick)
axs[1,1].plot(xs,ds['nexcite6'],label='c',linewidth=lthick)
axs[1,1].legend()
#axs[1,1].set_xscale('log')
axs[1,1].set_yscale('log')
axs[1,1].set_xlim([xmin-vs,xmax-vs])
axs[1,1].set_xlabel('$x_s$')
axs[1,1].set_ylabel('Level populations')


"""axs[1,1].plot(ds['xgrid']/ds['time'],ds['ro_p'],label='ro_p',color='b',linewidth=lthick)
axs[1,1].plot(ds['xgrid']/ds['time'],ds['ro_n'],label='ro_n',color='r',linewidth=lthick)
axs[1,1].plot(ds['xgrid']/ds['time'],ds['ro_n']+ds['ro_p'],label='ro_n+ro_p',color='g',linewidth=lthick)
axs[1,1].plot(dsm['xgrid']/dsm['time'],dsm['ro_p'],label='ro_p (MHD)',color='k',linewidth=lthick)
axs[1,1].legend()
axs[1,1].set_xscale('log')
#axs[1,1].set_yscale('log')
axs[1,1].set_xlim([0.01,2.0])
axs[1,1].set_xlabel('$x/t$')
axs[1,1].set_ylabel('Density')"""

plt.savefig(sname,dpi=300)
#plt.savefig('shocksub2_plot_midc.png',dpi=300)
#plt.savefig('test_plot.png',dpi=300)

T1=ds['pr_p'][-1]/ds['ro_p'][-1]*5.0/6.0
T=ds['pr_p']/ds['ro_p']*5.0/6.0
tpost=np.interp(xmin,ds['xgrid']/ds['time'],T*T0/T1)
tpre=np.interp(xmax,ds['xgrid']/ds['time'],T*T0/T1)
tmax=np.max(T[1000:-1]*T0/T1)
ropost=np.interp(xmin,ds['xgrid']/ds['time'],ds['ro_p']+ds['ro_n'])
ropre=np.interp(xmax,ds['xgrid']/ds['time'],ds['ro_p']+ds['ro_n'])

print(tpre,tpost,tmax,ropost/ropre)

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
