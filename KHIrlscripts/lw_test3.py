#Forward modelling test

#from lightweaver.fal import Falc82
from lightweaver.rh_atoms import H_6_atom, C_atom, O_atom, Si_atom, Al_atom, \
CaII_atom, Fe_atom, He_9_atom, MgII_atom, N_atom, Na_atom, S_atom
#from Si29_gkerr_cmo_updated import Si29
import lightweaver as lw
import matplotlib.pyplot as plt
import time
import numpy as np
import pipreadmods as PIPpy
from pipatmos import PIP1d 

#lw.benchmark()

#stop

def synth_lines(atmos, conserve, useNe, wave,activeatoms):
    '''
    Synthesise a spectral line for given atmosphere with different
    conditions.

    Parameters
    ----------
    atmos : lw.Atmosphere
        The atmospheric model in which to synthesise the line.
    conserve : bool
        Whether to start from LTE electron density and conserve charge, or
        simply use from the electron density present in the atomic model.
    useNe : bool
        Whether to use the electron density present in the model as the
        starting solution, or compute the LTE electron density.
    wave : np.ndarray
        Array of wavelengths over which to resynthesise the final line
        profile for muz=1.

    Returns
    -------
    ctx : lw.Context
        The Context object that was used to compute the equilibrium
        populations.
    Iwave : np.ndarray
        The intensity at muz=1 for each wavelength in `wave`.
    '''
    # Configure the atmospheric angular quadrature
    atmos.quadrature(5)
    # Configure the set of atomic models to use.
    aSet = lw.RadiativeSet([H_6_atom(), C_atom(), O_atom(), Si_atom(),
                            Al_atom(), CaII_atom(), Fe_atom(), He_9_atom(),
                            MgII_atom(), N_atom(), Na_atom(), S_atom()
                           ])
    # Set H and Ca to "active" i.e. NLTE, everything else participates as an
    # LTE background.
    aSet.set_active('H',activeatoms)
    # Compute the necessary wavelength dependent information (SpectrumConfiguration).
    spect = aSet.compute_wavelength_grid()

    # Either compute the equilibrium populations at the fixed electron density
    # provided in the model, or iterate an LTE electron density and compute the
    # corresponding equilibrium populations (SpeciesStateTable).
    if useNe:
        eqPops = aSet.compute_eq_pops(atmos)
    else:
        eqPops = aSet.iterate_lte_ne_eq_pops(atmos)

    # Configure the Context which holds the state of the simulation for the
    # backend, and provides the python interface to the backend.
    # Feel free to increase Nthreads to increase the number of threads the
    # program will use.
    ctx = lw.Context(atmos, spect, eqPops, conserveCharge=conserve, Nthreads=1)
    # Iterate the Context to convergence (using the iteration function now
    # provided by Lightweaver)
    lw.iterate_ctx_se(ctx)
    # Update the background populations based on the converged solution and
    # compute the final intensity for mu=1 on the provided wavelength grid.
    eqPops.update_lte_atoms_Hmin_pops(atmos)
    Iwave = ctx.compute_rays(wave, [atmos.muz[-1]], stokes=False)
    return ctx, Iwave

#fname='../KHIrldata/CD_800_1.2/'
#ds=PIPpy.pipread(fname,70)

#pipatmos.pip1d(ds,400)

#activeatoms='Ca'
#activeatoms='Si'
activeatoms='Mg'

#wave = np.linspace(853.9444, 854.9444, 1001) #For Ca II 8542
#wave = np.linspace(140.3, 140.4, 1001) #For Si IV 1403.77
wave = np.linspace(278.0, 286.0, 1001) #For Mg II

#atmosRef = PIP1d()
#ctxRef, IwaveRef = synth_lines(atmosRef, conserve=False, useNe=True, wave=wave,activeatoms=activeatoms)
atmosCons = PIP1d()
ctxCons, IwaveCons = synth_lines(atmosCons, conserve=True, useNe=False, wave=wave,activeatoms=activeatoms)
#atmosLte = PIP1d()
#ctx, IwaveLte = synth_lines(atmosLte, conserve=False, useNe=False, wave=wave,activeatoms=activeatoms)

#plt.plot(wave, IwaveRef, label='Reference FAL')
plt.plot(wave, IwaveCons, label='Reference Cons')
#plt.plot(wave, IwaveLte, label='Reference LTE n_e')
plt.legend()
plt.show()


#lw.atmosphere.ScaleType(0)

#lw.Atmosphere.make_1d(scale=lw.atmosphere.ScaleType(0),depthScale=ds['ygrid'],vlos=ds['vy_p'][1,:],vturb=0.0*ds['vy_p'][1,:],
#		temperature=temperature[:,1],ne=ds['ro_p'][:,1])
#lw.Atmosphere.make_2d()

#plt.contourf(np.log10(temperature),levels=101)
#cb.colorbar

#plt.show()
