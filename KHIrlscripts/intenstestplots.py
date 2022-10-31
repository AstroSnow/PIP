#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 15:08:10 2022

Following:
    https://chiantipy.readthedocs.io/en/latest/quick_start.html#g-n-t-function

@author: ben
"""
import ChiantiPy.core as ch
import numpy as np
import matplotlib.pyplot as plt
import pipreadmods#_old as pipreadmods
from scipy.interpolate import interp1d

#Load the spectral line
temperature = np.logspace(3.0,7.0,101)
si4 = ch.ion('si_4',temperature=temperature,eDensity=1.e+9,em=1.e+27)
mg2 = ch.ion('mg_2',temperature=temperature,eDensity=1.e+9,em=1.e+27)
#si4.popPlot()
#plt.show()

#Plot intensity vs temperature
si4.intensity()
#dist = np.abs(np.asarray(si4.Intensity['wvl']) - 1393.75)
dist = np.abs(np.asarray(si4.Intensity['wvl']) - 1402.77)
idx = np.argmin(dist)
print(' wvl = %10.3f '%(si4.Intensity['wvl'][idx]))
plt.loglog(temperature,si4.Intensity['intensity'][:,idx])
f = interp1d(temperature,si4.Intensity['intensity'][:,idx], kind='cubic',fill_value=0.0,bounds_error=False)

mg2.intensity()
#dist = np.abs(np.asarray(si4.Intensity['wvl']) - 1393.75)
dist = np.abs(np.asarray(mg2.Intensity['wvl']) - 2795.53)
idx = np.argmin(dist)
print(' wvl = %10.3f '%(mg2.Intensity['wvl'][idx]))
plt.loglog(temperature,mg2.Intensity['intensity'][:,idx])
f2 = interp1d(temperature,mg2.Intensity['intensity'][:,idx], kind='cubic',fill_value=0.0,bounds_error=False)

#Load in simulation data
filename = "/home/ben/Documents/KHI/PIPrl/KHIrldata/CD_800_1.2/"
#filename = '/home/ben/Documents/KHI/simdata/KHI_2D_r0/'
intensarr_si4=np.empty([71,812])
intensarr_mg2=np.empty([71,812])
tarr=np.empty([71])

for t in range(10,70):
    ds=pipreadmods.pipread(filename,tstep=t,vararrin=['pr_p','ro_p'])
    xg=ds['xgrid']
    yg=ds['ygrid']
    tarr[t]=ds['time']
    
    simtemp=1.0e6*ds['pr_p']/ds['ro_p']*5.0/3.0 #Dimensionalise the temperature
    emm_si4=simtemp*0.0 #create an empty array for the emmision
    emm_mg2=simtemp*0.0 #create an empty array for the emmision
    for i in range(0,np.size(xg)):
        for j in range(0,np.size(yg)):
            emm_si4[j,i]=ds['ro_p'][j,i]**2*f(simtemp[j,i])
            emm_mg2[j,i]=ds['ro_p'][j,i]**2*f2(simtemp[j,i])
    
    intens_siv=np.sum(emm_si4,axis=1)
    intensarr_si4[t,:]=intens_siv
    intens_mg2=np.sum(emm_mg2,axis=1)
    intensarr_mg2[t,:]=intens_mg2

fig,(ax,ax2)=plt.subplots(2,1,dpi=300)
ax.set_facecolor('k')
fig.set_size_inches(9.7,6.0)
#plt.xlabel('Time')
#plt.ylabel('Total mass below $T_{cut}$')
#line1, = ax.plot(yg,intens_siv)
cvels=np.linspace(-6,6,101)
cvels2=np.linspace(0,6,101)
cp=ax.contourf(tarr[10:70],yg,np.log10(intensarr_si4[10:70,:].T),levels=cvels2,cmap='magma',extend='both')
ax2.set_facecolor('k')
cp2=ax2.contourf(tarr[10:70],yg,np.log10(intensarr_mg2[10:70,:].T),levels=cvels,cmap='inferno',extend='both')
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.82, 0.55, 0.02, 0.4])
fig.colorbar(cp, cax=cbar_ax,label='log(Si IV intensity)')
cbar_ax = fig.add_axes([0.82, 0.08, 0.02, 0.4])
fig.colorbar(cp2, cax=cbar_ax,label='log(Mg IIk intensity)')

savename=('Figures/intensplots_rl_h.png')
plt.savefig(savename)