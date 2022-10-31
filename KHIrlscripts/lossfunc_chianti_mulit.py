#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 10:28:30 2022

@author: ben
"""
import ChiantiPy.core as ch
import numpy as np
import matplotlib.pyplot as plt
import pipreadmods_old as pipreadmods
import h5py

#minabund=1.e-4 #H, He, C, O, Ne
#minabund=2.e-5 #H, He, C, O, Ne, N, Mg, Si, S
minabund=1.e-6 #H, He, C, O, Ne, N, Mg, Si, S, Na, Al, Ar, Ca, Ni

#abundfile='sun_coronal_2012_schmelz_ext.abund'

nelements=101

temp = np.logspace(3.0,7.0,nelements)
rl = ch.radLoss(temp, 1.e+4, minAbund=minabund,verbose=1)

hf = h5py.File('lossfunc_default.h5', 'w')

hf.create_dataset('temperature', data=temp)
hf.create_dataset('rad_loss', data=rl.RadLoss['rate'])
hf.create_dataset('minAbund', data=minabund)
#hf.create_dataset('Abundances', data=abundfile)
hf.create_dataset('nelements', data=nelements)
hf.close()

#########################################################################

abundfile='sun_coronal_2012_schmelz_ext.abund'

nelements=101

temp = np.logspace(3.0,7.0,nelements)
rl2 = ch.radLoss(temp, 1.e+4, minAbund=minabund,abundance=abundfile,verbose=1)

hf = h5py.File('lossfunc_coronal_schmelz_ext.h5', 'w')

hf.create_dataset('temperature', data=temp)
hf.create_dataset('rad_loss', data=rl2.RadLoss['rate'])
hf.create_dataset('minAbund', data=minabund)
hf.create_dataset('Abundances', data=abundfile)
hf.create_dataset('nelements', data=nelements)
hf.close()


########################################################################

abundfile='sun_photospheric_2015_scott.abund'

nelements=101

temp = np.logspace(3.0,7.0,nelements)
rl3 = ch.radLoss(temp, 1.e+4, minAbund=minabund,abundance=abundfile,verbose=1)

hf = h5py.File('lossfunc_photo_scott.h5', 'w')

hf.create_dataset('temperature', data=temp)
hf.create_dataset('rad_loss', data=rl3.RadLoss['rate'])
hf.create_dataset('minAbund', data=minabund)
hf.create_dataset('Abundances', data=abundfile)
hf.create_dataset('nelements', data=nelements)
hf.close()

########################################################################

abundfile='sun_photospheric_2011_caffau.abund'

nelements=101

temp = np.logspace(3.0,7.0,nelements)
rl4 = ch.radLoss(temp, 1.e+4, minAbund=minabund,abundance=abundfile,verbose=1)

hf = h5py.File('lossfunc_photo_caffau.h5', 'w')

hf.create_dataset('temperature', data=temp)
hf.create_dataset('rad_loss', data=rl4.RadLoss['rate'])
hf.create_dataset('minAbund', data=minabund)
hf.create_dataset('Abundances', data=abundfile)
hf.create_dataset('nelements', data=nelements)
hf.close()


###################################################################

fig,ax=plt.subplots(1,1,dpi=300)
fig.set_size_inches(9.7,6.0)
plt.xlabel('log(T)')
plt.ylabel('Loss rate')
line1, = ax.plot(temp,rl.RadLoss['rate'],linewidth=3,label='Default')
ax.set_xscale('log')
ax.set_yscale('log')
line2, = ax.plot(temp,rl2.RadLoss['rate'],linewidth=3,label='Coronal_schemlz_ext')
line3, = ax.plot(temp,rl3.RadLoss['rate'],linewidth=3,label='photospheric_scott')
line4, = ax.plot(temp,rl4.RadLoss['rate'],linewidth=3,label='photospheric_caffau')
#line6, = ax.plot(tarrhd,tengarrhd, label='HD (time/10)')
ax.legend(handles=[line1,line2,line3,line4])

savename=('Figures/lossfun_chianti_multi_plot.png')
plt.savefig(savename)

