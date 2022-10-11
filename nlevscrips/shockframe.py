#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 08:47:02 2022
Alfven Mach numbers for the MHD (black) and plasma (blue) based on the bulk density. The red line shows the neutral sonic Mach number.
Requires simulation data

@author: ben
"""

import pipreadmods as PIPpy
import matplotlib.pyplot as plt
import numpy as np

#fname='../simdata/col_t10000/' ;sname='col_10000' #xi_n=0.94787198761607150
fname='../simdata/rad_thin_t10000/' ;sname='rad_thin_10000' #xi_n=0.94787198863788269
#fname='../simdata/rad_thick_t10000/' ;sname='rad_thick_10000' #xi_n=0.94787198863788269

ds=PIPpy.pipread(fname,10)

dsm=PIPpy.pipread('../simdata/MHD_t10000/',10) # Reference MHD simulation

#a=min(deriv(ds.vx_p),b)
a=np.gradient(ds['vx_p'])
b=np.argmin(a)

vs=ds['xgrid'][b]/ds['time']
xs=ds['xgrid']/ds['time']-vs


#am=min(deriv(dsm.vx_p),bm)
#vsm=dsm.x(bm)/dsm.t(0)
#xsm=dsm.x/dsm.t(0)-vsm
am=np.gradient(dsm['vx_p'])
bm=np.argmin(am)

vsm=dsm['xgrid'][bm]/dsm['time']
xsm=dsm['xgrid']/dsm['time']-vsm


#a=min(deriv(ds.vx_n),bn)                          
#vsn=ds.x(bn)/ds.t(0) 
an=np.gradient(ds['vx_n'])
bn=np.argmin(an)

vsn=ds['xgrid'][bn]/ds['time']



va=ds['bx']/np.sqrt(ds['ro_p']+ds['ro_n'])
vam=dsm['bx']/np.sqrt(dsm['ro_p'])

csn=np.sqrt(5.0/3.0*ds['pr_n']/ds['ro_n'])

fig, axs = plt.subplots(1,1)
fig.set_size_inches(9.7, 6)

lthick=2.5

axs.plot(xsm,abs(dsm['vx_p']-vsm)/vam,'k--',label='MHD',linewidth=lthick)
axs.plot(xs,abs(ds['vx_p']-vs)/va,'b',label='Alfven',linewidth=lthick)
axs.plot(xs,abs(ds['vx_n']-vsn)/csn,'r',label='Sonic',linewidth=lthick)

axs.set_xlabel('$x_s$')
axs.set_ylabel('Mach number')
axs.set_xlim([-0.015,0.01])
axs.legend()

plt.savefig(''.join(['shockframe_',sname]),dpi=300)