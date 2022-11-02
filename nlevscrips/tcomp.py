#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 14:53:51 2022

@author: ben
"""

import pipreadmods as PIPpy
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

#fname='../simdata/T_6220_long/Data/' ;xmin=0.0748;xmax=0.08 ;sname='shocksub2_plot_upc.png'; etime=30;T0=6220
#fname='../simdata/isca_midc_ltest/Data/' ;xmin=0.052;xmax=0.074 ;sname='shocksub2_plot_midc.png'; etime=40; T0=5030
#fname='../simdata/isca_lowc_long/Data/' ;xmin=0.04;xmax=0.1 ;sname='shocksub2_plot_lowc.png'; etime=21; T0=5180
#fname='../simdata/T_5300_n_7.5e16/' ;xmin=0.05;xmax=0.08 ;sname='shocksub2_plot_T5300.png'; etime=40;T0=5300
matplotlib.rcParams.update({'font.size': 18})

#ds=PIPpy.pipread(fname,etime)
dsm=PIPpy.pipread('../simdata/MHD_ref/',50)

nelements=7
tprearr=np.empty(nelements)
tpostarr=np.empty(nelements)
tmaxarr=np.empty(nelements)
rarr=np.empty(nelements)
xnarr=np.empty(nelements)
tref=np.empty(nelements)

for i in range(0,nelements):
    if i==0:
#        fname='../simdata/isca_lowc_long/Data/' ;xmin=0.04;xmax=0.1; etime=21; T0=5180
#    if i==1:
        fname='../simdata/isca_midc_ltest/Data/' ;xmin=0.052;xmax=0.074; etime=40; T0=5030
    if i==1:
        fname='../simdata/T_5300_n_7.5e16/' ;xmin=0.05;xmax=0.08; etime=40;T0=5300
    if i==2:
        fname='../simdata/T_5600_n_7.5e16/' ;xmin=0.05;xmax=0.068; etime=80;T0=5600
    if i==3:
        fname='../simdata/T_6000_n_7.5e16/' ;xmin=0.05;xmax=0.06; etime=40;T0=6000
    if i==4:
        fname='../simdata/T_6100_n_7.5e16/' ;xmin=0.06;xmax=0.064; etime=40;T0=6100
    if i==5:
        fname='../simdata/T_6220_long/Data/' ;xmin=0.0748;xmax=0.08; etime=30;T0=6220
    if i==6:
        fname='../simdata/T_6500_n_7.5e16/' ;xmin=0.141;xmax=0.144; etime=40;T0=6500
    
    ds=PIPpy.pipread(fname,etime)
    T1=ds['pr_p'][-1]/ds['ro_p'][-1]*5.0/6.0
    T=ds['pr_p']/ds['ro_p']*5.0/6.0
    tpost=np.interp(xmin,ds['xgrid']/ds['time'],T*T0/T1)
    tpre=np.interp(xmax,ds['xgrid']/ds['time'],T*T0/T1)
    tmax=np.max(T[100:-1]*T0/T1)
    ropost=np.interp(xmin,ds['xgrid']/ds['time'],ds['ro_p']+ds['ro_n'])
    ropre=np.interp(xmax,ds['xgrid']/ds['time'],ds['ro_p']+ds['ro_n'])
    xn=ds['ro_n']/(ds['ro_p']+ds['ro_n'])

    print(tpre,tpost,tmax,ropost/ropre,xn[-1])
    tprearr[i]=tpre
    tpostarr[i]=tpost
    tmaxarr[i]=tmax
    rarr[i]=ropost/ropre
    xnarr[i]=xn[-1]
    tref[i]=T0

dsm=PIPpy.pipread('../simdata/MHD_ref/',50)
T1m=dsm['pr_p'][-1]/dsm['ro_p'][-1]*5.0/3.0
Tm=dsm['pr_p']/dsm['ro_p']*5.0/3.0
tpost=np.interp(0.05,dsm['xgrid']/dsm['time'],Tm*T0/T1m)
tpre=np.interp(0.25,dsm['xgrid']/dsm['time'],Tm*T0/T1m)
tmax=np.max(Tm*T0/T1m)
ropost=np.interp(0.05,dsm['xgrid']/dsm['time'],dsm['ro_p'])
ropre=np.interp(0.25,dsm['xgrid']/dsm['time'],dsm['ro_p'])

tprem=tpre
tpostm=tpost
tmaxm=tmax
rm=ropost/ropre
xnm=0.0
trefm=T0
#xn=np.array([0.87,0.9986,0.9997,0.999992])
#xi=1.0-xn

#mhdtjump=27772.766819280467/5526.754752598247

#tpre=np.array([6158.830514202272,5237.838535379454,4894.209028189453,4927.655727575849])
#tpost=np.array([8394.621226583751,7862.661566436549,7936.2773700190355,9048.333527801136])
#tmax=np.array([24222.73214722533,16817.874630670856,17006.501840154626,17835.3943491788])

#comp=np.array([4.464929983051906,6.053326213145618,5.983905519969008,5.806633285098448])

fig, axs = plt.subplots(1, 1,dpi=300)
fig.set_size_inches(9.7, 6)
lthick=2.5

"""line1=plt.plot(np.log10(1.0-xnarr),tpostarr/tprearr,'-*r')
line2=plt.plot(np.log10(1.0-xnarr),tmaxarr/tprearr,'-*b')
line3=plt.plot(np.log10(1.0-xnarr),rarr,'-*g')"""

line1=plt.plot(tref,tpostarr/tprearr,'-*r',label='$T^d/T^u$',linewidth=lthick)
line2=plt.plot(tref,tmaxarr/tprearr,'-*b',label='$T_{max}/T^u$',linewidth=lthick)
line3=plt.plot(tref,rarr,'-*g',label='$\\rho^d/\\rho^u$',linewidth=lthick)
line4=plt.plot([5000,6500],[tpostm/tprem,tpostm/tprem],'r--',linewidth=lthick)#,label='$T^d/T^u$ (MHD)'
#axs=plt.plot([5000,6500],[tmaxm/tprem,tmaxm/tprem],'b--')
line5=plt.plot([5000,6500],[rm,rm],'g--',linewidth=lthick)#,label='$\\rho^d/\\rho^u$ (MHD)'
axs.legend()
#axs[1,1].set_xscale('log')
#axs[1,1].set_yscale('log')
#axs[1,1].set_xlim([xmin-vs,xmax-vs])
axs.set_xlabel('$Temperature [K]$')
axs.set_ylabel('$T_{jump}, \\rho_{jump}$')
plt.savefig('tcomp_plot',dpi=300)

#print(xi)