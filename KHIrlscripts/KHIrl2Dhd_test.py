#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 09:38:03 2023

@author: ben
"""
import numpy as np
import matplotlib.pyplot as plt
import pipreadmods


#filename = "/home/bs428/Documents/PIPrl/KHI2D_hr/"
filename = "/home/bs428/Documents/PIPrl/KHI2Dr1e2_hr/"

#tstep=[1,2,3,4,5,6,7]

fig,(ax)=plt.subplots(1,3,dpi=300)
fig.set_size_inches(9.7,6.0)

ii=0

cvelsT=np.linspace(3.5,6.01,101)
cvelsro=np.linspace(-0.05,2.1,101)
cvelsloss=np.logspace(-6,-2,101)

ds=pipreadmods.pipread(filename,tstep=12)

T=ds['pr_p']/ds['ro_p']*5.0/3.0*1.0e6
cp=ax[0].contourf(ds['xgrid'],ds['ygrid'],np.log10(T),levels=cvelsT,cmap='plasma')
cp2=ax[1].contourf(ds['xgrid'],ds['ygrid'],np.log10(ds['ro_p']),levels=cvelsro,cmap='gray')
cp3=ax[2].contourf(ds['xgrid'],ds['ygrid'],ds['edref_m'],levels=cvelsloss,cmap='twilight')

    
savename=('KHI2Dhr_test_r1e2.png')
plt.savefig(savename)