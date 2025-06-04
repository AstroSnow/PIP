import numpy as np
import astropy.constants as const
import astropy.units as u
import lightweaver as lw
from lightweaver.rh_atoms import H_9_atom
from scipy.integrate import quad

lambda_edges = [c.lambdaEdge for c in H_9_atom().continua]

# def integrand_numerator(nu, Trad):
#     return 1 / (np.exp(const.h.value * nu / (const.k_B.value * Trad)) - 1)

# def integrand_denominator(nu, Trad):
#     return 1 / nu * 1 / (np.exp(const.h.value * nu / (const.k_B.value * Trad)) - 1)

def integrand_numerator_x(x, Trad):
    return 1 / (np.exp(x) - 1)

def integrand_denominator_x(x, Trad):
    return const.h.value / const.k_B.value / Trad / x * 1 / (np.exp(x) - 1)

def wave_to_x(wave, Trad):
    nu = (wave << u.nm).to(u.Hz, equivalencies=u.spectral()).value
    return const.h.value * nu / const.k_B.value / Trad

Trad = 5777
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
for lambda_edge, energy in zip(lambda_edges, avg_energies):
    print(f"{lambda_edge:.2f} nm: {energy.value:.2f} eV")
