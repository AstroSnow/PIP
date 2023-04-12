#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 09:38:03 2023

@author: ben
"""
import numpy as np
import matplotlib.pyplot as plt
import pipreadmods


filename = "/home/bs428/Documents/KHI3D/KHI3Dr1e4/"

#
#pltvar='ro_p'
pltvar='loss'

fig = plt.figure(dpi=300)
ax = fig.add_subplot(111, projection='3d')
#fig.set_size_inches(9.7,6.0)

cvelsro=np.linspace(-0.05,2.1,101)

ds=pipreadmods.pipread(filename,tstep=42,vararrin=['ro_p'])

#meanro=np.mean(np.mean(ds['ro_p'],axis=2),axis=0)

if pltvar == 'ro_p':
	#Density
	lro=np.transpose(np.log10(ds['ro_p']),[0,2,1])
	g=np.where(np.logical_and(lro > 0.99, lro < 1.01))
	c=lro[g[0],g[1],g[2]]
	ccol='gray'

if pltvar == 'loss':
	#Losses
	loss=np.transpose(ds['edref_m'],[0,2,1])
	g=np.where(loss > 8.0e-5)
	c=loss[g[0],g[1],g[2]]
	ccol='plasma'

#ax.scatter(g[0],g[1],g[2],c=c,s=1, alpha=0.2,marker='s')
ax.scatter((g[0]/(512.0/2.0)-1.0)*50.0,g[1]/203.0-1.0,2.0*(g[2]-309.0)/412.0,c=c,s=1, alpha=0.2,marker='s',cmap=ccol)
ax.view_init(elev=15., azim=40)
#ax.plot_surface(g[0],g[1],g[2], alpha=0.1)

ax.set_xlim([-50, 50])
ax.set_ylim([-1, 1])
ax.set_zlim([-1.5, 0.5])
  
#savename=('KHI3D_volumeplot.png')
#plt.savefig(savename)