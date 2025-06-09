import numpy as np
import astropy.constants as const
import astropy.units as u
import lightweaver as lw
from lightweaver.rh_atoms import H_9_atom
from scipy.integrate import quad
import matplotlib.pyplot as plt
import h5py

Trad = 5777
T_elec=6000.0 #electron temperature

lambda_edges = [c.lambdaEdge for c in H_9_atom().continua]

# def integrand_numerator(nu, Trad):
#     return 1 / (np.exp(const.h.value * nu / (const.k_B.value * Trad)) - 1)

# def integrand_denominator(nu, Trad):
#     return 1 / nu * 1 / (np.exp(const.h.value * nu / (const.k_B.value * Trad)) - 1)

def integrand_numerator_x(x, Trad):
    return T_elec * const.k_B.value/const.h.value * np.exp(-x) / (1.0-np.exp(-x*T_elec/Trad))

def integrand_denominator_x(x, Trad):
    return x * np.exp(-x) / (1.0-np.exp(-x*T_elec/Trad))

def wave_to_x(wave, Trad):
    nu = (wave << u.nm).to(u.Hz, equivalencies=u.spectral()).value
    return const.h.value * nu / const.k_B.value / T_elec

def get_energies(T_elec):
    avg_freqs = []
    for lambda_edge in lambda_edges:
        num_int, num_int_err = quad(
            integrand_numerator_x,
            wave_to_x(lambda_edge, Trad),
            np.inf,
            args=(Trad,),
            epsabs=1e-25
        )
        denom_int, denom_int_err = quad(
            integrand_denominator_x,
            wave_to_x(lambda_edge, Trad),
            np.inf,
            args=(Trad,),
            epsabs=1e-25
        )

        if num_int_err / abs(num_int) > 1e-5 or denom_int_err / abs(denom_int) > 1e-5:
            print(f"Error may be high for lambda_edge {lambda_edge:.2f} nm: {num_int_err:.3e}, {denom_int_err:.3e}")

        avg_freqs.append(num_int / denom_int)

    avg_energies = [(const.h * (f << u.Hz)).to(u.eV) for f in avg_freqs]
    #for lambda_edge, energy in zip(lambda_edges, avg_energies):
    #    print(f"{lambda_edge:.2f} nm: {energy.value:.2f} eV")
    return avg_energies
    
n_elements=1001
T_elec_arr=np.logspace(3,6,n_elements)
cooling_energy=np.zeros_like(T_elec_arr)
cooling_energy1=np.zeros_like(T_elec_arr)
cooling_energy2=np.zeros_like(T_elec_arr)
cooling_energy3=np.zeros_like(T_elec_arr)
cooling_energy4=np.zeros_like(T_elec_arr)
cooling_energy5=np.zeros_like(T_elec_arr)
for i in range(0,np.size(T_elec_arr)):
    T_elec=T_elec_arr[i]
    avg_energies=get_energies(T_elec)
    print(T_elec,(avg_energies[0].value))
    cooling_energy[i]=avg_energies[0].value
    cooling_energy1[i]=avg_energies[1].value
    cooling_energy2[i]=avg_energies[2].value
    cooling_energy3[i]=avg_energies[3].value
    cooling_energy4[i]=avg_energies[4].value
    cooling_energy5[i]=avg_energies[5].value

#for i in range(1,n_elements):
#    dLt=np.log(T_elec_arr[i])-np.log(T_elec_arr[i-1])
#    print(dLt)

#Save the data    
f = h5py.File("ave_photon_energy_rec.hdf5", "w")
dset = f.create_dataset("n_elements",data=n_elements)
dset = f.create_dataset("T_elec",data=np.log(T_elec_arr))
dset = f.create_dataset("p-n0",data=cooling_energy)
dset = f.create_dataset("p-n1",data=cooling_energy1)
dset = f.create_dataset("p-n2",data=cooling_energy2)
dset = f.create_dataset("p-n3",data=cooling_energy3)
dset = f.create_dataset("p-n4",data=cooling_energy4)
dset = f.create_dataset("p-n5",data=cooling_energy5)
f.close()
#plt.plot(np.log10(T_elec_arr),np.log10(cooling_energy)) 
#plt.plot(np.log10(T_elec_arr),4.9-np.log10(T_elec_arr),'r')
#plt.show()

plt.loglog(T_elec_arr,cooling_energy,label='p-n0')
plt.loglog(T_elec_arr,cooling_energy1,label='p-n1')
plt.loglog(T_elec_arr,cooling_energy2,label='p-n2')
plt.loglog(T_elec_arr,cooling_energy3,label='p-n3')
plt.loglog(T_elec_arr,cooling_energy4,label='p-n4')
plt.loglog(T_elec_arr,cooling_energy5,label='p-n5')
plt.legend()
plt.show() 
