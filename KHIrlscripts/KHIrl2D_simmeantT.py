#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 09:38:03 2023

@author: ben
"""
import numpy as np
import matplotlib.pyplot as plt
import pipreadmods


filename = "/home/ben/Documents/KHI/PIPrl/KHIrldata/KHIrl2D1e2/"

fig,(ax)=plt.subplots(1,1,dpi=300)
fig.set_size_inches(9.7,6.0)

meanT=np.zeros(15)
time=np.zeros(15)

for i in range(0,15):

    ds=pipreadmods.pipread(filename,tstep=i)

    T=ds['pr_p']/ds['ro_p']*5.0/3.0*1.0e6
    time[i]=ds['time']
    meanT[i]=np.mean(T)

cp=ax.plot(time,np.log10(meanT))
plt.xlabel('time')
plt.ylabel('$log_{10}(\hat{T}$)')
    
savename=('../Figures2D/KHIrl2D_simmeanT.png')
plt.savefig(savename)