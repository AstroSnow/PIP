from typing import Callable

import numpy as np

import lightweaver.constants as Const
from lightweaver.atmosphere import Atmosphere, ScaleType

import pipreadmods as PIPpy
import matplotlib.pyplot as plt

fname='../KHIrldata/CD_800_1.2/'

ds=PIPpy.pipread(fname,70)

height=400

temperature=ds['pr_p']/ds['ro_p']*5.0/3.0*1.0e6 #Temperature in Kelvin


cmass = ds['ro_p'][height,:]*10**-5 #in KG/M^2?

temp = temperature[height,:] 

#plt.plot(np.log10(temp))
#plt.show()

ne = ds['ro_p'][height,:]*1.0e10#/(1.0e-30) #in m^-3

nh = np.outer([1.3575E+05, 6.8537E-02, 4.1543E-02, 5.2083E-02, 7.0588E-02, 1.0457E+10],ds['ro_p'][height,:]*1.0e10) #in m^-3

vel = ds['vx_p'][height,:]*8.0 #in km/s

vturb = ds['vy_p'][height,:]*1.0e1# np.zeros_like(ne) #in m/s 2-5km/s

#set upper and lower bc as zero radiation

PIP1d: Callable[[], Atmosphere] = lambda: Atmosphere.make_1d(ScaleType.ColumnMass, depthScale=cmass*Const.G_TO_KG / Const.CM_TO_M**2, temperature=np.copy(temp), ne=ne / Const.CM_TO_M**3, vlos=vel * Const.KM_TO_M, vturb=vturb * Const.KM_TO_M, hydrogenPops=nh / Const.CM_TO_M**3)
#, hydrogenPops=nh / Const.CM_TO_M**3)
