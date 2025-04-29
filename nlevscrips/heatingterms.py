#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 25 08:05:55 2022

@author: ben
"""
import pipreadmods as PIPpy
import matplotlib.pyplot as plt
import numpy as np

#fname='../simdata/col_t10000/' ;sname='col_10000' #xi_n=0.94787198761607150
#fname='../simdata/rad_thin_t10000/' ;sname='rad_thin_10000' #xi_n=0.94787198863788269
fname='../simdata/rad_thick_t10000/' ;sname='rad_thick_10000' #xi_n=0.94787198863788269

ds=PIPpy.pipread(fname,10)

tp=ds['pr_p']/ds['ro_p']*5.0/6.0
tn=ds['pr_n']/ds['ro_n']*5.0/3.0

beta=0.1

#######################################################################
#Collisional coupling
ac=np.sqrt(0.5*(tp+tn))/np.sqrt(beta/2.0*5.0/3.0)

#######################################################################
#Frictional heating
fricheat=ac*ds['ro_n']*ds['ro_p']*(np.sqrt(ds['vx_n']**2+ds['vy_n']**2)-np.sqrt(ds['vx_p']**2+ds['vy_p']**2))**2

#######################################################################
#Thermal damping
tdamp=3.0/2.0*ac*ds['ro_p']*ds['ro_n']*(ds['pr_n']/ds['ro_n']-ds['pr_p']/ds['ro_p']/2.0)

tdampp=tdamp
np.where(tdampp > 0, tdampp, 0)

tdampn=-tdamp
np.where(tdampn > 0, tdampn, 0)

#######################################################################
#Heating/cooling
hc=ds['aheat']-ds['ion_loss']

hcp=hc
np.where(hcp > 0, hcp, 0)

hcn=-hc
np.where(hcn > 0, hcn, 0)

######################################################################
#plots
fig, axs = plt.subplots(2,2)
fig.set_size_inches(9.7, 6)


#axs.plot(ds['xgrid'],ac,label='ac')
axs[0,0].plot(ds['xgrid']/ds['time'],fricheat,label='$F_{heat}$',color='k')
axs[0,0].plot(ds['xgrid']/ds['time'],tdampp,label='$T_{damp} > 0$',color='b')
axs[0,0].plot(ds['xgrid']/ds['time'],tdampn,label='$T_{damp} > 0$',color='r')
axs[0,0].set_yscale('log')
#axs[0].set_xscale('log')
axs[0,0].set_ylim([10**-20,10**-2])
axs[0,0].set_xlim([0.01,0.05])
axs[0,0].legend(loc='lower right')

#axs[1,0].plot(ds['xgrid']/ds['time'],tn-tp,label='T_n-T_p')
axs[1,0].plot(ds['xgrid']/ds['time'],ds['ro_n'],label='$\\rho _{n} $',color='r')
axs[1,0].plot(ds['xgrid']/ds['time'],ds['ro_p'],label='$\\rho _{p} $',color='b')
#axs[1].set_xscale('log')
axs[1,0].set_xlim([0.01,0.05])
axs[1,0].legend()

axs[0,1].plot(ds['xgrid']/ds['time'],hcp,label='$A_{heat} > \phi_{ion}$',color='g')
axs[0,1].plot(ds['xgrid']/ds['time'],hcn,label='$A_{heat} < \phi_{ion}$',color='orange')
axs[0,1].set_yscale('log')
#axs[2].set_xscale('log')
axs[0,1].set_ylim([10**-7,10**-2])
axs[0,1].set_xlim([0.01,0.05])
axs[0,1].legend()

axs[1,1].plot(ds['xgrid']/ds['time'],ds['vx_n'],label='$v_{xn}$',color='r')
axs[1,1].plot(ds['xgrid']/ds['time'],ds['vx_p'],label='$v_{xp}$',color='b')
#axs[1,1].set_yscale('log')
#axs[2].set_xscale('log')
#axs[1,1].set_ylim([10**-7,10**-2])
axs[1,1].set_xlim([0.01,0.05])
axs[1,1].legend()

plt.savefig(''.join(['heating_terms',sname]),dpi=300)