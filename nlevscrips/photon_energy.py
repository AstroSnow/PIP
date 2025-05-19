#Calculate the mean photon energy
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

c=299792458.0 #m/s
h=6.62607015e-34 #J s
heV=4.1357e-15 # eV s
kB=1.380649e-23 #J/K
#kB=8.6173303e-5 #eV/K
Trad=5000.0
Te=8000.00

#Waveloength of transition
#l0=c/91.175e-9
l0=91.2e-9
l1=410.0e-9
#l0=91.175e-9
#print(l0,1.0/l0)

def photon_ionisation(Trad, Te,nu):
    res=1.0/nu*(1.0+1.0/(np.exp(h*nu/kB/Trad)-1.0))*np.exp(-h*nu/kB/Te)
    #print(np.exp(-h*nu/kB/Te))
    return(res)

def photon_recombination(Trad,nu):
    res=1.0/nu*1.0/(np.exp(h*nu/kB/Trad)-1.0)
    return(res)

def plank(nu):
    E=8.0*np.pi*heV*nu**3/c**3* 1.0/(np.exp(h*nu/kB/Trad)-1.0)
    return E

def plank_l(l):
    E=8.0*np.pi*heV*l**-3*1.0/(np.exp(h*c/l/kB/Trad)-1.0)
    return E
#nu=1.0/np.logspace(-3,5,1000)
#l=np.linspace(0.1*l0,100*l0,10000)

#print(h/l/kB/Trad)
#Spectral energy
#E=1.0/(np.exp(h*l/kB/Trad)-1.0)*l**3

#l=np.linspace(0.1,100,10000)*1.0e-15
#E=1.0/(np.exp(h/l/kB/Trad)-1.0)/l**3

#nu=np.linspace(l0,100000000*l0,10000)/c
l=np.linspace(1.0e-8,1.0e-5,10000)
nu=c/l
#E=8.0*np.pi*heV*nu**3/c**3* 1.0/(np.exp(h*nu/kB/Trad)-1.0)
#E=E*6.242e18 # convert joules to eV
#l=nu/c*1.0e-6 #in micrometres

#l=np.linspace(l0*0.001,100*l0,10000)
#E=1.0/(np.exp(h*l/kB/Trad)-1.0)*l**3

#plt.loglog(nu,1.0/(np.exp(h/nu/kB/Trad)-1.0))
#plt.plot(l,plank(nu),color='k')
plt.plot(l,plank_l(l),color='k')
plt.axvline(x=l0,color='b')
plt.axvline(x=l1,color='r')
#plt.plot([l1,l1],[0,6000])
#plt.loglog(1.0/nu,photon_recombination(Trad,nu)*4.13558e-15)
#plt.loglog(1.0/nu,photon_ionisation(Trad,Te,nu)*4.13558e-15)
plt.show()

print(integrate.quad(plank,nu[-1],nu[0]))
print(integrate.quad(plank_l,l0,l[-1]))

"""
#########################################################################
Notes
#########################################################################
Collisional Ionisation
#########################################################################
Let E be the energy of the Blackbody field
Let F be the integrand from sollum: F=1/nu 1/(exp(h nu / kB T) -1)
Then the mean energy the ionising photon is then
M= (int(E(nu) F(nu) d(nu) )/int(F(nu) d(nu)
 I.e., the energy of the field, multiplied by the probability, divided by the probability.
CHECK THAT E(nu0)=13.6 for ground state. i.e., that the minimum ionising photon has at least the ionising energy

The excess heating is then phi_(R,I)=Gamma_{R,I} (M-phi_0)/\hat{\phi}
where Gamma_{R,I} is the photoionisation rate, M is defined above, phi_0 is the energy of the transition, \hat{\phi} is the normalisation factor in the code

#########################################################################
Collisional Recombination
#########################################################################
Let E be the energy of the maxwellian distribution of electrons at the given temperature (Maxwell-Boltzman?)
Let F be the integrand from sollum: F=1/nu (1+1/(exp(h nu / kB T) -1))exp(-h nu/kB T_e)
Then the mean energy the recombining electron is then
M_R= (int(E(nu) F(nu) d(nu) )/int(F(nu) d(nu)
The cooling term is then phi_(R,R)=Gamma_{R,R} (M_R-phi_0)/\hat{\phi}
NOT SURE ABOUT THE phi_0 TERM.....

Questions: Collisional recombination includes both spontanious and induced terms. For the spontanious part, the photon energy loss is the energy from the electron field minus the energy required to be stored for the transition? What about induced?

"""
