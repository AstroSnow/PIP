#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 11:06:14 2022

@author: ben
"""
import numpy as np
import matplotlib.pyplot as plt
import pipreadmods


#filename = "/home/ben/Documents/KHI/PIPrl/KHIrldata/isca_3d_test/"
filename = "/home/ben/Documents/KHI/PIPrl/KHIrldata/isca_3d_test_2/"
#filename = "/home/ben/Documents/KHI/PIPrl/KHIrldata/isca_3d_test_zpert/Data/"
#filename = "/home/ben/Documents/KHI/PIPrl/KHIrldata/CD_800_1.2/"

ds=pipreadmods.pipread(filename,tstep=60,vararrin=['pr_p','ro_p'])

fig,(ax)=plt.subplots(1,1,dpi=300)
ax.set_facecolor('k')
fig.set_size_inches(9.7,6.0)
#plt.xlabel('Time')
#plt.ylabel('Total mass below $T_{cut}$')
#line1, = ax.plot(yg,intens_siv)
cvels=np.linspace(-6,6,101)
#cvels2=np.linspace(0,6,101)
T=ds['pr_p']/ds['ro_p']*5.0/3.0*1.0e6
#cp=ax.contourf(ds['xgrid'],ds['ygrid'],np.log10(T[175,:,:]),levels=101)
cp=ax.contourf(ds['xgrid'],ds['zgrid'],np.log10(T[:,280,:]),levels=101)
#cp=ax.contourf(ds['xgrid'],ds['ygrid'],np.log10(T[50,:,:]),levels=101)
#cp=ax.contourf(ds['xgrid'],ds['zgrid'],np.log10(T[:,60,:]),levels=101)
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.82, 0.1, 0.02, 0.8])
fig.colorbar(cp, cax=cbar_ax,label='log(T)')

savename=('Figures/plottest3d.png')
plt.savefig(savename)
