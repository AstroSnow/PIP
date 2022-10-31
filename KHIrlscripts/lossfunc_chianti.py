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
minabund=2.e-5 #H, He, C, O, Ne, N, Mg, Si, S
#minabund=1.e-6 #H, He, C, O, Ne, N, Mg, Si, S, Na, Al, Ar, Ca, Ni

abundfile='sun_coronal_2012_schmelz_ext.abund'

nelements=101

temp = np.logspace(3.0,7.0,nelements)
rl = ch.radLoss(temp, 1.e+4, minAbund=minabund,abundance=abundfile,verbose=1)

rl.radLossPlot()

hf = h5py.File('lossfunc.h5', 'w')

hf.create_dataset('temperature', data=temp)
hf.create_dataset('rad_loss', data=rl.RadLoss['rate'])
hf.create_dataset('minAbund', data=minabund)
hf.create_dataset('Abundances', data=abundfile)
hf.create_dataset('nelements', data=nelements)

hf.close()