#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 09:38:03 2023

@author: ben
"""
import numpy as np
import matplotlib.pyplot as plt
import pipreadmods


filename = "/home/bs428/Documents/PIPrl/KHI2D_hr/"

tstep=[1,2,3,4,5,6,7]

fig,(ax)=plt.subplots(2,7,dpi=300)
fig.set_size_inches(21.0,6.0)

ii=0

cvelsT=np.linspace(3.5,6.01,101)
cvelsro=np.linspace(-0.05,2.1,101)

for i in tstep:

    ds=pipreadmods.pipread(filename,tstep=i)

    T=ds['pr_p']/ds['ro_p']*5.0/3.0*1.0e6
    cp=ax[0,ii].contourf(ds['xgrid'],ds['ygrid'],np.log10(T),levels=cvelsT,cmap='plasma')
    cp2=ax[1,ii].contourf(ds['xgrid'],ds['ygrid'],np.log10(ds['ro_p']),levels=cvelsro,cmap='gray')
    
    ii=ii+1
    
savename=('KHI2D_hr_context.png')
plt.savefig(savename)