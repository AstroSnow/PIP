#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 14:53:51 2022

@author: ben
"""

import pipreadmods as PIPpy
import matplotlib.pyplot as plt
import numpy as np

xn=np.array([0.87,0.9997,0.999992])
xi=1.0-xn

mhdtjump=27772.766819280467/5526.754752598247

tpre=np.array([6158.830514202272,4894.209028189453,4927.655727575849])
tpost=np.array([8394.621226583751,7936.2773700190355,9048.333527801136])
tmax=np.array([24222.73214722533,17006.501840154626,17835.3943491788])

comp=np.array([4.464929983051906,5.983905519969008,5.806633285098448])

fig, axs = plt.subplots(1, 1,dpi=300)
fig.set_size_inches(9.7, 6)
lthick=2.5

axs=plt.plot(np.log10(xi),tpost/tpre,'-*')
axs=plt.plot(np.log10(xi),tmax/tpre,'-*')
axs=plt.plot(np.log10(xi),comp,'-*')

print(xi)