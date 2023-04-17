#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 08:37:37 2023

@author: bs428
"""

import pipreadmods as PIPpy
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.animation as manimation

n=11
fname = "/home/bs428/Documents/PIPrl/KHI2Dr1e2_hr/"
ds=PIPpy.pipread(fname,0)


# Define the meta data for the movie
FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='KHI2Drl_hr_mov', artist='Ben Snow',
                comment='KHI simulation with radiative losses')
writer = FFMpegWriter(fps=2, metadata=metadata)


fig,(axs)=plt.subplots(1,3,dpi=100)
fig.set_size_inches(12, 4)

cvelsT=np.linspace(3.5,6.01,101)
cvelsro=np.linspace(-0.05,2.1,101)

divider = make_axes_locatable(axs[0])
cax = divider.append_axes('right', size='5%', pad=0.05)
axs[0].set_aspect('equal')
T=ds['pr_p']/ds['ro_p']*5.0/3.0*1.0e6
cf=axs[0].contourf(ds['xgrid'],ds['ygrid'],np.log10(T),levels=cvelsT,cmap='plasma')
fig.colorbar(cf, cax=cax, orientation='vertical')

divider = make_axes_locatable(axs[1])
cax2 = divider.append_axes('right', size='5%', pad=0.05)
axs[1].set_aspect('equal')
cf=axs[1].contourf(ds['xgrid'],ds['ygrid'],np.log10(ds['ro_p']),levels=cvelsro,cmap='gray')
fig.colorbar(cf, cax=cax2, orientation='vertical')

divider = make_axes_locatable(axs[2])
cax3 = divider.append_axes('right', size='5%', pad=0.05)
axs[2].set_aspect('equal')
cp3=axs[2].contourf(ds['xgrid'],ds['ygrid'],ds['edref_m'])
fig.colorbar(cf, cax=cax3, orientation='vertical')

# Update the frames for the movie
with writer.saving(fig, "KHI2Drl_hr_mov_r1e2.mp4", 100):
    for i in range(n):
        ds=PIPpy.pipread(fname,i)
        
        #fig, axs = plt.subplots(1, 1,dpi=300)
        #fig.set_size_inches(9.7, 6)
        #divider = make_axes_locatable(axs)
        #cax = divider.append_axes('right', size='5%', pad=0.05)
        #axs.set_aspect('equal')

        T=ds['pr_p']/ds['ro_p']*5.0/3.0*1.0e6
        cf=axs[0].contourf(ds['xgrid'],ds['ygrid'],np.log10(T),levels=cvelsT,cmap='plasma')        
        
        cf=axs[1].contourf(ds['xgrid'],ds['ygrid'],np.log10(ds['ro_p']),levels=cvelsro,cmap='gray')
        fig.colorbar(cf, cax=cax2, orientation='vertical')
        
        cp3=axs[2].contourf(ds['xgrid'],ds['ygrid'],ds['edref_m'])
        fig.colorbar(cf, cax=cax3, orientation='vertical')
    
        writer.grab_frame()

#plt.show()