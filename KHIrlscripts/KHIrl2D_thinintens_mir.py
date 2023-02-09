#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 09:59:51 2023

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
filename = "../KHIrldata/KHIrl2D1e2/"
#filename = '/home/ben/Documents/KHI/simdata/KHI_2D_r0/'
tsamps=15
#gsamps=823#
gsamps=412
intensarr_si4=np.empty([tsamps,gsamps])
intensarr_mg2=np.empty([tsamps,gsamps])
tarr=np.empty([tsamps])

for t in range(0,15):
    ds=pipreadmods.pipread(filename,tstep=t,vararrin=['pr_p','ro_p'])
    xg=ds['xgrid']
    yg=ds['ygrid']
    tarr[t]=ds['time']
    
    simtemp=1.0e6*ds['pr_p']/ds['ro_p']*5.0/3.0 #Dimensionalise the temperature
    
    A=ds['ro_p']
    B= np.flipud(A)
    ro= np.concatenate((A[1:],B), axis=0)
    A=simtemp
    B= np.flipud(A)
    simtemp= np.concatenate((A[1:],B), axis=0)
    
    emm_si4=simtemp*0.0 #create an empty array for the emmision
    emm_mg2=simtemp*0.0 #create an empty array for the emmision
    for i in range(0,np.size(xg)):
        for j in range(0,np.size(yg)*2-1):
            emm_si4[j,i]=ro[j,i]**2*max(f(simtemp[j,i]),0)
            emm_mg2[j,i]=ro[j,i]**2*max(f2(simtemp[j,i]),0)
    
    intens_siv=np.sum(emm_si4,axis=0)
    intensarr_si4[t,:]=intens_siv
    intens_mg2=np.sum(emm_mg2,axis=0)
    intensarr_mg2[t,:]=intens_mg2
    

fig,(ax)=plt.subplots(1,1,dpi=300)
#ax.set_facecolor('k')
fig.set_size_inches(9.7,6.0)
#plt.xlabel('Time')
#plt.ylabel('Total mass below $T_{cut}$')
#line1, = ax.plot(yg,intens_siv)
si4_mean=np.mean(intensarr_si4.T,axis=0)
mg2_mean=np.mean(intensarr_mg2.T,axis=0)
p1=ax.plot(tarr,si4_mean/si4_mean[0])
p2=ax.plot(tarr,mg2_mean/mg2_mean[0])
ax.set_yscale('log')
#cvels=np.linspace(-6,6,101)
#cvels2=np.linspace(0,6,101)
#cp=ax.contourf(np.log10(intensarr_si4[:,:].T),levels=cvels2,cmap='magma',extend='both')
#ax2.set_facecolor('k')
#cp2=ax2.contourf(np.log10(intensarr_mg2[:,:].T),levels=cvels,cmap='inferno',extend='both')
#fig.subplots_adjust(right=0.8)
#cbar_ax = fig.add_axes([0.82, 0.55, 0.02, 0.4])
#fig.colorbar(cp, cax=cbar_ax,label='log(Si IV intensity)')
#cbar_ax = fig.add_axes([0.82, 0.08, 0.02, 0.4])
#fig.colorbar(cp2, cax=cbar_ax,label='log(Mg IIk intensity)')

#savename=('Figures/intensplots_rl_h.png')
#plt.savefig(savename)