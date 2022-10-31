#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 10:17:21 2022

@author: ben
"""
import numpy as np
import matplotlib.pyplot as plt
import pipreadmods


filename = "/home/ben/Documents/KHI/PIPrl/KHIrldata/isca_3d_test/"
#filename = "/home/ben/Documents/KHI/PIPrl/KHIrldata/CD_800_1.2/"

meanT=np.empty([15,412])
meanL=np.empty([15,412])
meanP=np.empty([15,412])
tarr=np.empty(15)

for i in range(0,15):
    ds=pipreadmods.pipread(filename,tstep=i*10)
    T=ds['pr_p']/ds['ro_p']*5.0/3.0*1.0e6
    meanT[i,:]=(np.mean(np.mean(T,axis=0),axis=1))
    meanL[i,:]=(np.mean(np.mean(ds['edref_m'],axis=0),axis=1))
    meanP[i,:]=(np.mean(np.mean(ds['pr_p']/(5.0/3.0-1.0),axis=0),axis=1))
    tarr[i]=ds['time']
    
fig,(ax)=plt.subplots(1,1,dpi=300)
#ax.set_facecolor('k')
fig.set_size_inches(9.7,6.0)
#plt.xlabel('Time')
#plt.ylabel('Total mass below $T_{cut}$')
#line1, = ax.plot(yg,intens_siv)
#cvels=np.linspace(-6,6,101)
#cvels2=np.linspace(0,6,101)
#T=ds['pr_p']/ds['ro_p']*5.0/3.0*1.0e6
#cp=ax.contourf(ds['xgrid'],ds['ygrid'],np.log10(T[175,:,:]),levels=101)
#cp=ax.contourf(ds['ygrid'],tarr,np.log10(meanT),levels=101)
#fig.subplots_adjust(right=0.8)
#cbar_ax = fig.add_axes([0.82, 0.1, 0.02, 0.8])
#fig.colorbar(cp, cax=cbar_ax,label='log(T)')
plt.plot(tarr,meanP[:,0])
plt.xlabel('Time')
plt.ylabel('mean thermal energy at y=0')
plt.plot([0,100],[meanP[0,0],meanP[0,0]-100.0*0.003])

savename=('Figures/loss3dtest.png')
plt.savefig(savename)
