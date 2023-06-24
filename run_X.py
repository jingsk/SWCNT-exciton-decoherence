import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from SWCNT import SWCNT

" CONSTANTS "

# Boltzmann constant (J/K)
K_b = 1.380649e-23
# planck's constant
h_bar = 1.054571817e-34 #Js

#time to evaluate Chi
t_array = np.linspace(0, 3, num=100) * 1e-12 #s




#unitless vibrational coupling
def gamma(q):
    G = D_s * q / np.sqrt(2 * rho * L * h_bar * v_s * q)
    F = np.exp(-(q ** 2 * sigma ** 2)/4)
    w = v_s * q
    return (G * F / w)
#phonon occupation number
def n(q):
    return (np.exp((h_bar * v_s * q)/(K_b * T)) - 1) ** -1

#the function to be integrated
def X_integrand(q,t):
    return gamma(q)**2 *(1+2*n(q)) * (1-np.cos(v_s * q *t))

def X0_integrand(q,t):
    return gamma(q)**2 * (1-np.cos(v_s * q *t))

def XT_integrand(q,t):
    return gamma(q)**2 *(2*n(q)) * (1-np.cos(v_s * q *t))

def norm(integ):
    exponents = -L / (1 * np.pi) * np.array([quad(integ, 0, 1*np.pi/L, args=(t))[0] for t in t_array])
    norm = np.exp(exponents)
    return norm

if __name__ == "__main__":
    n,m = 6,4
    swcnt = SWCNT(n,m)
    L = swcnt.L
    sigma = swcnt.L
    rho = swcnt.rho
    v_s = swcnt.v_s
    D_s = swcnt.D_s
    #time to evaluate Chi
    t_array = np.linspace(0, 3, num=100) * 1e-12 #s
    #simulation temperature in K
    T = 4 #K
    norms = [norm(X_integrand),norm(XT_integrand), norm(X0_integrand)]

    t_array2=t_array*1e12
    fig, ax = plt.subplots(figsize=(5,4))
    ax.plot(t_array2,norms0_1[0], t_array2,norms0_1[1], t_array2,norms0_1[2])
    #ax.axis('off')
    axs.set_xlabel('Time (ps)')
    plt.setp(axs, yscale='log', ylim=[0.08,1.1],xlim=[-0.01,2.5])
    ax.xlabel('Time (ps)')
    ax.legend(['|X|','|X_T|','|X_0|'],fontsize=7.5)
    plt.tight_layout()
    plt.savefig('X.png',dpi=900)
