#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 14:58:48 2022

@author: ben
"""
#Saha equilibrium test for 6 level hydrogen

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
#from mpl_toolkits.axes_grid1 import host_subplot
#import mpl_toolkits.axisartist as AA
from mpl_toolkits.axisartist.parasite_axes import HostAxes, ParasiteAxes

matplotlib.rcParams.update({'font.size': 18})

#A few VALC datapoint to test (T0 in K, nelec in cm^-3
harr= [ -25.0,   0.0,  50.0, 100.0, 150.0,
        250.0, 350.0, 450.0, 515.0, 555.0,
        605.0, 655.0, 705.0, 755.0, 855.0,
        905.0, 980.0,1065.0,1180.0,1280.0,
       1380.0,1515.0,1605.0,1785.0,1925.0,
       1990.0,2016.0,2050.0,2070.0,2080.0,
       2090.0,2104.0,2107.0,2109.0,2113.0,
       2115.0,2120.0,2129.0,2160.0,2200.0,
       2230.0,2255.0]

T0arr=[6910.0,6420.0,5840.0,5455.0,5180.0,
       4780.0,4465.0,4220.0,4170.0,4230.0,
       4420.0,4730.0,5030.0,5280.0,5650.0,
       5755.0,5925.0,6040.0,6150.0,6220.0,
       6280.0,6370.0,6440.0,6630.0,6940.0,
       7160.0,7360.0,7660.0,7940.0,8180.0,
       8440.0,9500.0,10700.0,12300.0,18500.0,
       21000.0,22500.0,23000.0,23500.0,24000.0,
       24200.0,24500.0]


nelecarr=[1.547e14,6.433e13,2.122e13,1.066e13,6.476e12,
          2.674e12,1.110e12,4.516e11,2.495e11,1.733e11,
          1.112e11,8.085e10,7.664e10,8.838e10,1.064e11,
          1.049e11,1.041e11,9.349e10,8.108e10,7.486e10,
          7.600e10,6.456e10,6.005e10,4.771e10,4.028e10,
          3.858e10,3.811e10,3.792e10,3.783e10,3.780e10,
          3.799e10,3.705e10,3.535e10,3.306e10,2.620e10,
          1.943e10,1.881e10,1.812e10,1.120e10,2.009e10,
          1.943e10,1.881e10] 



narr=np.zeros((6,np.size(harr)))
nuin=np.zeros(np.size(harr))
nuni=np.zeros(np.size(harr))
            
for j in range(0, 41):
    T0=T0arr[j]
    nelec=nelecarr[j]
    
    #Convert nelec to m^-3
    nelec=nelec*1.0e6
    
    
    #https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwji-JeZ_6bzAhUDmFwKHaCtBaYQFnoECAQQAQ&url=https%3A%2F%2Fwww.mtholyoke.edu%2Fcourses%2Fjlevine%2Fast228%2FHW%2Fboltzsahaprobsolutions.pdf&usg=AOvVaw3j0oply73jMUvsKx7vNnsg
    #T0=8000.0
    #nelec=20.0/kb/T0
    #nelec=0.012*0.7/kb/T0
    
    
    #constants
    kb=1.38064852e-23 #Boltzman constant in m^2 kg s^-2 K^-1
    pi=3.14 #pi
    mh=1.67e-27 #hydrogen mass in kg
    a0=5.29e-11 #Bohr radius in m
    me=9.10938356e-31 #Electron mass in kg
    h=6.62607004e-34 #Planck's constant in m2 kg s^-1
    
    
    n=np.zeros(7) #population density of different hydrogen types
    n[0]=0 #blank so I can port it to fortran easier
    n[1]=0 #ground state
    n[2]=0 #1st level excitation
    n[3]=0 #2nd level excitation
    n[4]=0 #3rd level excitation
    n[5]=0 #4th level excitation
    n[6]=1 #Ionised hydrogen
    
    g=np.zeros(7) #statistical weight of different hydrogen types
    g[0]=0 #blank so I can port it to fortran easier
    g[1]=2.0*1.0**2 #ground state
    g[2]=2.0*2.0**2 #1st level excitation
    g[3]=2.0*3.0**2 #2nd level excitation
    g[4]=2.0*4.0**2 #3rd level excitation
    g[5]=2.0*5.0**2 #4th level excitation
    g[6]=1.0 #Ionised hydrogen
    
    E=[0.0,13.6*1.6022e-19,3.4*1.6022e-19,1.51*1.6022e-19,0.85*1.6022e-19,0.54*1.6022e-19,0.0] #in eV 
    #E=E*1.6022e-19 #Convert to joules (to be dimensionally correct)

    
    n[6]=nelec
    n[1]=(2.0/nelec*g[6]/g[1]*(2.0*pi*me*kb*T0/h/h)**(3.0/2.0)*np.exp(-E[1]/kb/T0))
    n[2]=(2.0/nelec*g[6]/g[2]*(2.0*pi*me*kb*T0/h/h)**(3.0/2.0)*np.exp(-E[2]/kb/T0))
    n[3]=(2.0/nelec*g[6]/g[3]*(2.0*pi*me*kb*T0/h/h)**(3.0/2.0)*np.exp(-E[3]/kb/T0))
    n[4]=(2.0/nelec*g[6]/g[4]*(2.0*pi*me*kb*T0/h/h)**(3.0/2.0)*np.exp(-E[4]/kb/T0))
    n[5]=(2.0/nelec*g[6]/g[5]*(2.0*pi*me*kb*T0/h/h)**(3.0/2.0)*np.exp(-E[5]/kb/T0))
    
    #n(1)=nelec*g(1)/g(6)/2.0*(2.0*!pi*mehat*kbhat*T0/hhat/hhat*1.0e14)^(-3.0/2.0)*exp((E(6)-E(1))/kb/T0)
    
    n[1]=1.0/n[1]
    n[2]=1.0/n[2]
    n[3]=1.0/n[3]
    n[4]=1.0/n[4]
    n[5]=1.0/n[5]
    n[6]=n[6]/nelec
    
    ntot=n[1]+n[2]+n[3]+n[4]+n[5]+n[6]
    
    n=n/ntot
    
    print(n[1:7])
    
    narr[:,j]=n[1:7]
    
    #Coupling frequencies
    ni=nelec
    xii=(narr[5,j]/(narr[5,j]+narr[0,j]+narr[1,j]+narr[2,j]+narr[3,j]+narr[4,j])) #ion fraction
    nn=ni/xii-ni
    sigin=1.16e-18 #ion neutral cross section
    massin=mh/2.0
    
    nuin[j]=nn*sigin*np.sqrt(8*kb*T0/np.pi/massin)
    nuni[j]=ni*sigin*np.sqrt(8*kb*T0/np.pi/massin)

fig = plt.figure(dpi=300)
fig.set_size_inches(9.7, 6)
lthick=2.5

host = fig.add_axes([0.15, 0.1, 0.75, 0.8], axes_class=HostAxes)
par1 = ParasiteAxes(host, sharex=host)
host.parasites.append(par1)

host.axis["right"].set_visible(False)

par1.axis["right"].set_visible(True)
par1.axis["right"].major_ticklabels.set_visible(True)
par1.axis["right"].label.set_visible(True)

p1, = host.plot(harr,narr[5,:]/(narr[5,:]+narr[0,:]+narr[1,:]+narr[2,:]+narr[3,:]+narr[4,:]), label="Ion fraction",linewidth=lthick,color='k')
p2, = par1.plot(harr,T0arr, label="Temperature",linewidth=lthick,color='r')

p3, =host.plot([1380,1380],[0,2],'--') #mid-chromosphere test
p3, =host.plot([705,705],[0,2],'--')
p3, =host.plot([150,150],[0,2],'--')

host.set_yscale('log')
par1.set_yscale('log')
host.set_ylabel("Ion fraction")
host.set_xlabel('Height [km]')
par1.set_ylabel("Temperature")

host.set_ylim([1.0e-8,1.0])
#par1.set_ylim([3.0e3,4.0e4])

par1.axis["right"].label.set_color(p2.get_color())

#plt.savefig('saha2_plot.png',dpi=300)


####################################################

fig,axs = plt.subplots(1, 1,dpi=300)
fig.set_size_inches(9.7, 6)
lthick=2.5

axs.plot(harr,nuin,'k')
axs.plot(harr,nuni,'r')
axs.set_yscale('log')

"""fig, axs = plt.subplots(1, 1,dpi=300)
fig.set_size_inches(9.7, 6)
lthick=2.5
axs.plot(harr,narr[5,:]/(narr[5,:]+narr[0,:]+narr[1,:]+narr[2,:]+narr[3,:]+narr[4,:]),linewidth=lthick)
axs.set_yscale('log')
axs.set_xlabel('Height [km]')
axs.set_ylabel('Ion fraction')"""
