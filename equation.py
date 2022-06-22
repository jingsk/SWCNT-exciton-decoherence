import numpy as np
import warnings
import matplotlib.pyplot as plt


warnings.filterwarnings('ignore')

" CONSTANTS "

# move constant to a different files

# Boltzmann constant (J/K)
k_b = 1.380649e-23
# h bar (J-s)
h_bar = 1.054571817e-34
# temperature (K)
T = 4
# rho (kg/m)
p = 1.67e-15
# vs (m/s)
v = 19900
# deformation potential (J)
Ds = 1.9228e-18
# length of nanotube (m)
L = 40e-10
# confinement length (m)
sigma = 13.5e-10
# energy of the excited state (Hz)
omega_constant = 1.27 * 1.60218e-19 / h_bar
# Unknown constant from directly proportional relation (still working)
C = 1

# Chi - linear susceptibility, polarizability units are (C * m^2 / V*2)
# change Temperature to time
# make sure the omega is in Hertz so the exponential is unit less

" FUNCTIONS "

# wave vector
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

# temperature_dependent right side equation
temperature_dependent = 1j * np.exp((np.absolute(gamma) ** 2) * (- n * np.absolute(np.exp(- 1j * omega_s * t) - 1)) ** 2)
# extract imaginary part
temperature_dependent_imaginary = [ele.imag for ele in temperature_dependent]

plt.scatter(t, temperature_dependent_imaginary)
plt.ylabel('X(t) - Temperature Dependent')
plt.xlabel('Time(s)')
plt.show()

# temperature_independent right side equation
temperature_independent = 1j * np.exp(((np.absolute(gamma)) ** 2) * np.exp(- 1j * omega_s * t) - 1)
# extract imaginary part
temperature_independent_imaginary = [ele.imag for ele in temperature_independent]

plt.scatter(t, temperature_independent_imaginary)
plt.ylabel('X(t) - Temperature Independent')
plt.xlabel('Time(s)')
plt.show()

# phonon shifted transition frequency
big_omega = omega_constant * h_bar - np.absolute(gamma) ** 2 * omega_s

# linear susceptibility
linear_susceptibility = - 1j * np.exp(- 1j * big_omega * t) * C * temperature_independent * temperature_independent
# extract imaginary part
linear_susceptibility_imaginary = [ele.imag for ele in linear_susceptibility]

plt.scatter(t, linear_susceptibility_imaginary)
plt.ylabel('X(t)')
plt.xlabel('Time(s)')
plt.show()
