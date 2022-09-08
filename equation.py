import numpy as np
import warnings
import matplotlib.pyplot as plt
from scipy.integrate import quad

warnings.filterwarnings('ignore')

" CONSTANTS "

# Boltzmann constant (J/K)
K_b = 1.380649e-23
# planck's constant
h_bar = 1.054571817e-34 #Js
# length density (kg/m)
rho = 1.67e-15
# sound velocity (m/s)
v_s = 19900
# deformation potential (J)
D_s = 2.243047191e-18
# length of nanotube (m)
L = 18.83e-10
# confinement length (m)
sigma = 3e-9 #m
# energy of the excited state (Hz)
omega_constant = 1.41 * 1.60218e-19 / h_bar # Hz

#time to evaluate Chi
t_array = np.linspace(0, 3, num=1000) * 1e-12 #s

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
def integrand(q,t):
    return gamma(q)**2 *(1+2*n(q)) * (1-np.cos(v_s * q *t))

#high Temp
T = 6 #K
#the quad function integrate over q at time t
#the integral is integrated serially at t from t_array
#[0] after quad gives the value (see scipy.integrate doc)
exponents = -L / (1 * np.pi) * np.array([quad(integrand, 0, 1*np.pi/L, args=(t))[0] for t in t_array])
norm = np.exp(exponents)

#low Temp
T = 0.01 #K
exponents = -L / (1 * np.pi) * np.array([quad(integrand, 0, 1*np.pi/L, args=(t))[0] for t in t_array])
norm2 = np.exp(exponents)

#plotting
plt.plot(t_array,norm,t_array,norm2)
plt.yscale('log')
plt.title('|X(t)| vs. Time')
plt.ylabel('|X(t)|')
plt.xlabel('Time (ps)')
plt.legend(['6 K','0.1 K'])
plt.ylim([0.08,1.1])
plt.show()


# LONG TIME
t_array = np.linspace(0, 3, num=1000) * 1e-12 #s
T = 6

a = (D_s ** 2 * K_b * T) / (2 * h_bar ** 2 * v_s ** 3 * rho)
y = np.exp(- a * t_array)

plt.plot(t_array,y)
plt.ylim([0.08,1.1])
plt.yscale('log')
plt.title('Long time')
plt.ylabel('y = exp(- (Ds^2 Kb T)/(2 * hbar ** 2 * vs ** 3 * rho) t')
plt.xlabel('Time (ps)')
plt.show()

# CHIRALITY
T = 6

rho_low = 1.67e-15
rho_med = 3e-15
rho_high = 6e-15

a = (D_s ** 2 * K_b * T) / (2 * h_bar ** 2 * v_s ** 3 * rho_low)
y1 = np.exp(- a * t_array)

a = (D_s ** 2 * K_b * T) / (2 * h_bar ** 2 * v_s ** 3 * rho_med)
y2 = np.exp(- a * t_array)

a = (D_s ** 2 * K_b * T) / (2 * h_bar ** 2 * v_s ** 3 * rho_high)
y3 = np.exp(- a * t_array)

plt.plot(t_array, y1, t_array, y2, t_array,y3)
plt.ylim([0.08,1.1])
plt.yscale('log')
plt.title('|X(t)| (long time) vs Time (ps)')
plt.ylabel('|X(t)| (long time) at long time')
plt.xlabel('Time (ps)')
plt.legend(['(6,4)', '(12,8)', '(24,16)'])
plt.show()



