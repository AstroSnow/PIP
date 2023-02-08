#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 09:11:31 2023

@author: ben
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pipreadmods_old as pipreadmods
import h5py

#filename = "../lossfunc.h5"
filename='lossfunc_photo_scott.h5'

with h5py.File(filename, "r") as f:
    # Print all root level object names (aka keys) 
    # these can be group or dataset names 
    print("Keys: %s" % f.keys())
    # get first object name/key; may or may NOT be a group
    a_group_key = list(f.keys())[0]

    # get the object type for a_group_key: usually group or dataset
    print(type(f[a_group_key])) 

    # If a_group_key is a group name, 
    # this gets the object names in the group and returns as a list
    temperature = list(f['temperature'])
    loss=list(f['rad_loss'])

    # If a_group_key is a dataset name, 
    # this gets the dataset values and returns as a list
    #data = list(f[a_group_key])
    # preferred methods to get dataset values:
    #ds_obj = f[a_group_key]      # returns as a h5py dataset object
    #ds_arr = f[a_group_key][()]  # returns as a numpy array


loss2=np.zeros(np.size(loss))
loss3=np.zeros(np.size(loss))
for i in range(0,np.size(temperature)):
    loss2[i]=loss[i]
    loss3[i]=loss[i]
    if np.log10(temperature[i]) >= 5.5:
        loss2[i]=loss[i]*(1.0-np.tanh((np.log10(temperature[i])-5.5)/0.1)**2)
        loss3[i]=loss[i]*(1.0-np.tanh((np.log10(temperature[i])-5.5)/0.4)**2)

matplotlib.rcParams.update({'font.size': 20})
fig,ax=plt.subplots(1,1,dpi=300)
fig.set_size_inches(9.7,6.0)
plt.xlabel('log(T)')
plt.ylabel('Normalised loss')
line1, = ax.plot(np.log10(temperature),loss/np.max(loss),linewidth=3,label='Scott')
line2, = ax.plot(np.log10(temperature),loss2/np.max(loss2),linewidth=3,label='Modified')
#line3, = ax.plot(np.log10(temperature),loss3/np.max(loss3),linewidth=3,label='Modified (c=0.4)')
#line6, = ax.plot(tarrhd,tengarrhd, label='HD (time/10)')
ax.legend(handles=[line1,line2])

savename=('../Figures2D/KHIrl2D_lossprofile.png')
plt.savefig(savename)
