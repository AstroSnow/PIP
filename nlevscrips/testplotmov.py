#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 12:00:10 2022

@author: snow
"""

import pipreadmods as PIPpy
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.animation as manimation

n=160
fname='/home/snow/Documents/KHInlev/KHInlevdata/Data/' #;xmin=0.0748;xmax=0.08 ;sname='shocksub2_plot_upc.png'; etime=30;T0=6220
ds=PIPpy.pipread(fname,0)


# Define the meta data for the movie
FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='Movie Test', artist='Ben Snow',
                comment='KHI simulation with 6level H')
writer = FFMpegWriter(fps=2, metadata=metadata)


fig, axs = plt.subplots(2, 2,dpi=300)
fig.set_size_inches(9.7, 9.7)

divider = make_axes_locatable(axs[0,0])
cax = divider.append_axes('right', size='5%', pad=0.05)
axs[0,0].set_aspect('equal')
cf=axs[0,0].contourf(ds['ro_p'],levels=100,vmin=0.9,vmax=10.0)
fig.colorbar(cf, cax=cax, orientation='vertical')

divider = make_axes_locatable(axs[0,1])
cax2 = divider.append_axes('right', size='5%', pad=0.05)
axs[0,1].set_aspect('equal')
cf=axs[0,1].contourf(ds['ro_n'],levels=100,vmin=165,vmax=225)
fig.colorbar(cf, cax=cax2, orientation='vertical')

divider = make_axes_locatable(axs[1,0])
cax3 = divider.append_axes('right', size='5%', pad=0.05)
axs[1,0].set_aspect('equal')
cf=axs[1,0].contourf(ds['vx_p']-ds['vx_n'],levels=100,vmin=-0.003,vmax=0.003)
fig.colorbar(cf, cax=cax3, orientation='vertical')

divider = make_axes_locatable(axs[1,1])
cax4 = divider.append_axes('right', size='5%', pad=0.05)
axs[1,1].set_aspect('equal')
cf=axs[1,1].contourf(ds['vy_p']-ds['vy_n'],levels=100,vmin=-0.003,vmax=0.003)
fig.colorbar(cf, cax=cax4, orientation='vertical')

# Update the frames for the movie
with writer.saving(fig, "/home/snow/Documents/KHInlev/nlevscrips/KHnlev_data_noconv.mp4", 100):
    for i in range(n):
        ds=PIPpy.pipread(fname,i)
        
        #fig, axs = plt.subplots(1, 1,dpi=300)
        #fig.set_size_inches(9.7, 6)
        #divider = make_axes_locatable(axs)
        #cax = divider.append_axes('right', size='5%', pad=0.05)
        #axs.set_aspect('equal')

        cf=axs[0,0].contourf(ds['ro_p'],levels=100,vmin=0.9,vmax=10.0)
        fig.colorbar(cf, cax=cax, orientation='vertical')
        
        cf=axs[0,1].contourf(ds['ro_n'],levels=100,vmin=165,vmax=225)
        fig.colorbar(cf, cax=cax2, orientation='vertical')
        
        cf=axs[1,0].contourf(ds['vx_p']-ds['vx_n'],levels=100,vmin=-0.003,vmax=0.003)
        fig.colorbar(cf, cax=cax3, orientation='vertical')

        cf=axs[1,1].contourf(ds['vy_p']-ds['vy_n'],levels=100,vmin=-0.003,vmax=0.003)
        fig.colorbar(cf, cax=cax4, orientation='vertical')
        
        writer.grab_frame()

#plt.show()