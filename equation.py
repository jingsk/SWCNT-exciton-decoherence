import numpy as np
import math
import matplotlib.pyplot as plt


" CONSTANTS "

# Boltzmann constant (m2*kg*s-2*K-1)
k = 1.380649 * 10 ** -23
# h bar (J-s)
h = 1.054571817 * 10 ** -34
# temperature (K)
T = 4
# rho (kg/m)
p = 1.67 * 10 ** -15
# vs (km/s)
v = 19.9
# deformation potential (eV)
Ds = 1.9228 * 10 ** -18
# length of nanotube (m)
L = 10 ** -10
# confinement length (m)
sigma = 3.3 * 10 ** -10
# omega from the bare energy of the Quantum Dots state
omega_constant = 2


" FUNCTIONS "

# wave vector
q = np.linspace(0.1, 0.5, num=100)

# linear dispersion
omega = v * q

# deformation potential couplings
G = Ds * q / np.sqrt(2 * p * L * h * omega)

# form factor
F = np.exp(-(q ** 2 * sigma)/4)

# exciton coupling matrix elements
g = G * F

# linear dispersion (default values)
w = v * q

# dimensionless coupling strength
gamma = g/w

# phonon occupation number
n = (np.exp((h * w)/(k * T)) - 1) ** -1

# time_dependent right side equation
time_dependent = 1j * np.exp(((np.absolute(gamma)) ** 2) * (- n * np.absolute(np.exp(- 1j * omega * T) - 1)) ** 2)

# time_independent right side equation
time_independent = 1j * np.exp(((np.absolute(gamma)) ** 2) * np.exp(- 1j * omega * T) - 1)

# phonon shifted transition frequency
big_omega = omega_constant - np.absolute(gamma) ** 2 * omega

# linear susceptibility
X = - np.exp(- big_omega * T) * time_independent * time_independent

print(time_dependent)
# Here a sample on how to read a file
# sample_data = pd.read_csv('sample_data.csv')

