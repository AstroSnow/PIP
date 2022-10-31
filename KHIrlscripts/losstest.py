#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 10:17:21 2022

@author: ben
"""
import numpy as np
import matplotlib.pyplot as plt
import pipreadmods
from scipy.integrate import simps


filename = "/home/ben/Documents/KHI/PIPrl/KHIrldata/scott_loss_0.01/"
filename2 = "/home/ben/Documents/KHI/PIPrl/KHIrldata/scott_loss_0.001/"
filename3 = "/home/ben/Documents/KHI/PIPrl/KHIrldata/mod_loss_0.01/"

meanT=np.empty([15,412])
tarr=np.empty(15)
teng1=np.empty(15)

meanT2=np.empty([15,412])
tarr2=np.empty(15)
teng2=np.empty(15)

meanT3=np.empty([15,412])
tarr3=np.empty(15)
teng3=np.empty(15)

margin=6

for i in range(0,15):
    ds=pipreadmods.pipread(filename,tstep=i)
    T=ds['pr_p']/ds['ro_p']*5.0/3.0*1.0e6
    meanT[i,:]=np.mean(T,axis=1)
    tarr[i]=ds['time']
    teng=ds['pr_p']/(5.0/3.0-1.0)
    xg=ds['xgrid']
    yg=ds['ygrid']
    nx=len(xg)
    ny=len(yg)
    inteng=simps(simps(teng[margin:ny-margin-1,margin:nx-margin-1],xg[margin:nx-margin-1]),yg[margin:ny-margin-1])
    teng1[i]=inteng
    edr1=ds['edref_m']
    
    ds=pipreadmods.pipread(filename2,tstep=i)
    T2=ds['pr_p']/ds['ro_p']*5.0/3.0*1.0e6
    meanT2[i,:]=np.mean(T2,axis=1)
    tarr2[i]=ds['time']
    teng=ds['pr_p']/(5.0/3.0-1.0)
    xg=ds['xgrid']
    yg=ds['ygrid']
    nx=len(xg)
    ny=len(yg)
    inteng=simps(simps(teng[margin:ny-margin-1,margin:nx-margin-1],xg[margin:nx-margin-1]),yg[margin:ny-margin-1])
    teng2[i]=inteng
    edr2=ds['edref_m']
    
    ds=pipreadmods.pipread(filename3,tstep=i)
    T3=ds['pr_p']/ds['ro_p']*5.0/3.0*1.0e6
    meanT3[i,:]=np.mean(T3,axis=1)
    tarr3[i]=ds['time']
    teng=ds['pr_p']/(5.0/3.0-1.0)
    xg=ds['xgrid']
    yg=ds['ygrid']
    nx=len(xg)
    ny=len(yg)
    inteng=simps(simps(teng[margin:ny-margin-1,margin:nx-margin-1],xg[margin:nx-margin-1]),yg[margin:ny-margin-1])
    teng3[i]=inteng
    edr3=ds['edref_m']
    
cvels=np.linspace(3.5,6.1,101)
fig,((ax,ax2,ax5),(ax3,ax4,ax6))=plt.subplots(2,3,dpi=300)
#ax.set_facecolor('k')
fig.set_size_inches(9.7,6.0)
#plt.xlabel('Time')
#plt.ylabel('Total mass below $T_{cut}$')
#line1, = ax.plot(yg,intens_siv)
#cvels=np.linspace(-6,6,101)
#cvels2=np.linspace(0,6,101)
#T=ds['pr_p']/ds['ro_p']*5.0/3.0*1.0e6
cp=ax.contourf(ds['xgrid'],ds['ygrid'],np.log10(T),levels=cvels)
cp2=ax2.contourf(ds['xgrid'],ds['ygrid'],np.log10(T2),levels=cvels)
cp2=ax5.contourf(ds['xgrid'],ds['ygrid'],np.log10(T3),levels=cvels)
#cp=ax.contourf(ds['ygrid'],tarr,np.log10(meanT),levels=101)
#fig.subplots_adjust(right=0.8)
#cbar_ax = fig.add_axes([0.82, 0.1, 0.02, 0.8])
#fig.colorbar(cp, cax=cbar_ax,label='log(T)')
pl3=ax3.plot(tarr,meanT[:,0])
pl4=ax3.plot(tarr2,meanT2[:,0])
pl4=ax3.plot(tarr3,meanT3[:,0])
ax3.set_xlabel('Time')
ax3.set_ylabel('Mean coronal temperature')

pl5=ax4.plot(tarr,teng1)
pl5=ax4.plot(tarr2,teng2)
pl5=ax4.plot(tarr3,teng3)
ax4.set_xlabel('Time')
ax4.set_ylabel('Mean thermal energy')
#plt.ylabel('mean thermal energy at y=0')
#plt.plot([0,100],[meanP[0,0],meanP[0,0]-100.0*0.003])

pl5=ax6.plot((T.flatten()),edr1.flatten(),'.')
pl5=ax6.plot(T2.flatten(),edr2.flatten(),'.')
pl5=ax6.plot(T3.flatten(),edr3.flatten(),'.')
ax6.set_xlabel('T')
ax6.set_ylabel('Loss')
ax6.set_xscale('log')

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.82, 0.1, 0.02, 0.8])
fig.colorbar(cp, cax=cbar_ax,label='log(T)')

savename=('Figures/losstest.png')
plt.savefig(savename)


