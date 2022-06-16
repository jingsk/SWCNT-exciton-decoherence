import numpy as np
import math
import matplotlib.pyplot as plt


" CONSTANTS "

# move constant to a different files

# Boltzmann constant (m2*kg*s-2*K-1)
k_b = 1.380649e-23
# h bar (J-s)
h_bar = 1.054571817e-34
# temperature (K)
T = 4
# rho (kg/m)
p = 1.67e-15
# vs (km/s)
v = 19.9
# deformation potential (J)
Ds = 1.9228e-18
# length of nanotube (m)
L = 40e-15
# confinement length (m)
sigma = 13.5e-10
# energy of the excited state (Hz)
omega_constant = 1.27 * 1.60218e-19 / h_bar

# Chi - linear susceptibility, polarizability units are (C * m^2 / V*2)
# change Temperature to time
# make sure the omega is in Hertz so the exponential is unit less

" FUNCTIONS "

# wave vector
# q is a inverse distance (1/m)** Angstrom inverse, .1, .2 Angstrom
# print in a file
q = np.linspace(0.1, 0.5, num=100) * 1/L

# time (s)
t = np.linspace(0, 3, num=100) * 1e-12
# linear dispersion
omega_s = v * q

# deformation potential couplings
G = Ds * q / np.sqrt(2 * p * L * h_bar * omega_s)

# form factor
F = np.exp(-(q ** 2 * sigma ** 2)/4)

# exciton coupling matrix elements
g = G * F

# linear dispersion (default values)
w = v * q

# dimensionless coupling strength
gamma = g / w

# phonon occupation number
n = (np.exp((h_bar * w)/(k_b * T)) - 1) ** -1

# time_dependent right side equation
time_dependent = 1j * np.exp(((np.absolute(gamma)) ** 2) * (- n * np.absolute(np.exp(- 1j * omega_s * t) - 1)) ** 2)

# time_independent right side equation
time_independent = 1j * np.exp(((np.absolute(gamma)) ** 2) * np.exp(- 1j * omega_s * t) - 1)

# phonon shifted transition frequency
big_omega = omega_constant * h_bar - np.absolute(gamma) ** 2 * omega_s

# linear susceptibility
X = - 1j * np.exp(- 1j * big_omega * t) * time_independent * time_independent

print(q)
print(gamma)
print(time_dependent)
