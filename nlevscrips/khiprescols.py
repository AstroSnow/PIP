#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 15:19:42 2023

@author: bs428
"""
#Saha equilibrium test for 6 level hydrogen

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
#from mpl_toolkits.axes_grid1 import host_subplot
#import mpl_toolkits.axisartist as AA
from mpl_toolkits.axisartist.parasite_axes import HostAxes, ParasiteAxes
from scipy import interpolate
from operator import truediv

matplotlib.rcParams.update({'font.size': 18})

T0ref=5500.0
ro0ref=7.5e10
rojump=10.0
T0arr=[]
rojarr=[]
nelecarr=[ro0ref]
#nelecarr[1:-1]=ro0ref*rojump

#narr=np.zeros((6,2))
nuin=[]
nuni=[]
pnarr=[]
pparr=[]
ptarr=[]

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
       
j=0

################################
#Initial solution
T0arr=[T0ref]
T0=T0arr[0]
nelec=ro0ref

#Convert nelec to m^-3
nelec=nelec*1.0e6


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

#narr[:,j]=n[1:7]
narr0=[n[1]]
narr1=[n[2]]
narr2=[n[3]]
narr3=[n[4]]
narr4=[n[5]]
narr5=[n[6]]

#Coupling frequencies
ni=nelec
xii=(narr5[j]/(narr5[j]+narr0[j]+narr1[j]+narr2[j]+narr3[j]+narr4[j])) #ion fraction
xin=1.0-xii
nn=ni/xii-ni
sigin=1.16e-18 #ion neutral cross section
massin=mh/2.0

nuin=[nn*sigin*np.sqrt(8*kb*T0/np.pi/massin)]
nuni=[ni*sigin*np.sqrt(8*kb*T0/np.pi/massin)]

pn=nn*T0*3.0/5.0#xin/(xin+2.0*xii)
pp=ni*T0*6.0/5.0#2.0*xii/(xin+2.0*xii)
pt=pn+pp
pnarr=[pn]
pparr=[pp]
ptarr=[pt]

loopit=1

while loopit == 1:
    T0=T0+1.0
    j=j+1
    T0arr.append(T0)
    nelec=ro0ref*rojump
    
    #Convert nelec to m^-3
    nelec=nelec*1.0e6

    
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
    
    #narr[:,j]=n[1:7]
    narr0.append(n[1])
    narr1.append(n[2])
    narr2.append(n[3])
    narr3.append(n[4])
    narr4.append(n[5])
    narr5.append(n[6])
    
    #Coupling frequencies
    ni=nelec
    xii=(narr5[j]/(narr5[j]+narr0[j]+narr1[j]+narr2[j]+narr3[j]+narr4[j])) #ion fraction
    xin=1.0-xii
    nn=ni/xii-ni
    sigin=1.16e-18 #ion neutral cross section
    massin=mh/2.0
    
    nuin.append(nn*sigin*np.sqrt(8*kb*T0/np.pi/massin))
    nuni.append(ni*sigin*np.sqrt(8*kb*T0/np.pi/massin))

    pn=nn*T0*3.0/5.0#xin/(xin+2.0*xii)
    pp=ni*T0*6.0/5.0#2.0*xii/(xin+2.0*xii)
    pt=pn+pp
    pnarr.append(pn)
    pparr.append(pp)
    ptarr.append(pt)
    print(ptarr[j],ptarr[0],T0)
    
    if ptarr[j] < 0.5*ptarr[0]:
        loopit=0
    

#fig = plt.figure(dpi=300)
#fig.set_size_inches(9.7, 6)
#lthick=2.5


f=interpolate.interp1d(ptarr[1:-1],T0arr[1:-1],kind='cubic')
f2=interpolate.interp1d(T0arr[1:-1],ptarr[1:-1],kind='cubic')
res = list(map(truediv, narr5[1:-1], (narr5[1:-1]+narr0[1:-1]+narr1[1:-1]+narr2[1:-1]+narr3[1:-1]+narr4[1:-1])))
f3=interpolate.interp1d(T0arr[1:-1],res,kind='cubic')

T0sol=f(ptarr[0])
print(T0sol,f2(T0sol),ptarr[0],f3(T0sol))

#plt.savefig('saha2_plot.png',dpi=300)


