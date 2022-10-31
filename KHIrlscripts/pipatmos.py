from typing import Callable

import numpy as np

import lightweaver.constants as Const
from lightweaver.atmosphere import Atmosphere, ScaleType

import pipreadmods as PIPpy
import matplotlib.pyplot as plt

#fname='../KHIrldata/CD_800_1.2/'
fname = "/home/ben/Documents/KHI/PIPrl/KHIrldata/isca_3d_test_2/"
ds=PIPpy.pipread(fname,60)

xloc=200
zlog=200
L0=1.0e6 #Arbitrary guess for length scale as 1Mm

temperature=ds['pr_p']/ds['ro_p']*5.0/3.0*1.0e6 #Temperature in Kelvin


rho = ds['ro_p'][zloc,:,xloc]*10**-5 #in KG/M^2?
rho=np.append(rho,np.flip(rho))

temp = temperature[zloc,:,xloc] 
temp=np.append(temp,np.flip(temp))

#height coordinates?
depthscale=np.linspace(0,4,num=np.size(rho))*L0

#plt.plot(np.log10(temp))
#plt.show()

ne = ds['ro_p'][zloc,:,xloc]*1.0e10#/(1.0e-30) #in m^-3

#nh = np.outer([1.3575E+05, 6.8537E-02, 4.1543E-02, 5.2083E-02, 7.0588E-02, 1.0457E+10],ds['ro_p'][zloc,:,xloc]*1.0e10) #in m^-3
nHTot = rho / (lw.DefaultAtomicAbundance.massPerH * lw.Amu)


vel = ds['vz_p'][zloc,:,xloc]*8.0 #in km/s

vturb = np.zeros_like(ne)+2.0 #in m/s 2-5km/s #ds['vy_p'][height,:]*1.0e1# 

#set upper and lower bc as zero radiation

PIP1d: Callable[[], Atmosphere] = lambda: Atmosphere.make_1d(ScaleType.ColumnMass, depthScale=depthscale, temperature=np.copy(temp), ne=ne, vlos=vel * Const.KM_TO_M, vturb=vturb * Const.KM_TO_M, nHTot=nHTot,lowerBc: Optional[BoundaryCondition] = lw.ZeroRadiation(), upperBc: Optional[BoundaryCondition] = lw.ZeroRadiation())
#, hydrogenPops=nh / Const.CM_TO_M**3)
