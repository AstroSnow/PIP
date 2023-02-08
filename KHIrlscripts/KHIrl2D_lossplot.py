#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 09:08:21 2023

@author: ben
"""
import numpy as np
import matplotlib.pyplot as plt
import pipreadmods


#filename = "/home/ben/Documents/KHI/PIPrl/KHIrldata/isca_3d_test/"
filename = "/home/ben/Documents/KHI/PIPrl/KHIrldata/KHIrl2D1e2/"
#filename = "/home/ben/Documents/KHI/PIPrl/KHIrldata/isca_3d_test_zpert/Data/"
#filename = "/home/ben/Documents/KHI/PIPrl/KHIrldata/CD_800_1.2/"

ds=pipreadmods.pipread(filename,tstep=14)

fig,(ax)=plt.subplots(1,1,dpi=300)
ax.set_facecolor('k')
fig.set_size_inches(9.7,6.0)
T=ds['pr_p']/ds['ro_p']*5.0/3.0*1.0e6
cp=ax.plot(np.log10(T),ds['edref_m'],'.')
