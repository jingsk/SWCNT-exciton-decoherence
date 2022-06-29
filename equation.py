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
T = 1
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
omega_constant = 0.2 * 1.60218e-19 / h_bar
# Unknown constant from directly proportional relation
C = 1


" FUNCTIONS "


# wave vector
q = np.linspace(0.001, 0.01, num=100) * 1/L

# time (s)
t = np.linspace(0, 3, num=10000) * 1e-12
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
n_s = (np.exp((h_bar * w)/(k_b * T)) - 1) ** -1


# temperature_dependent
temperature_dependent = 1j * np.exp(np.absolute(gamma) ** 2 * - n_s
                                    @ (np.absolute(np.exp(- 1j * np.outer(omega_s, t)) - 1)) ** 2)

plt.scatter(t, temperature_dependent.imag)
plt.xlim(0, 2e-12)
plt.yscale('log')
plt.title('X_t vs Time(s)')
plt.ylabel('X_t')
plt.xlabel('Time(s)')
plt.show()

# temperature_independent
temperature_independent = 1j * np.exp(np.absolute(gamma) ** 2
                                      @ (np.exp(- 1j * np.outer(omega_s, t)) - 1))

plt.scatter(t, temperature_independent.imag)
plt.xlim(0, 2e-12)
plt.yscale('log')
plt.title('X_o vs Time (s)')
plt.ylabel('X_o')
plt.xlabel('Time(s)')
plt.show()

# phonon shifted transition frequency
big_omega = omega_constant * h_bar - np.absolute(gamma) ** 2 @ omega_s

# linear susceptibility
linear_susceptibility = C * - 1j * np.exp(- 1j * big_omega * t) \
                        * temperature_independent * temperature_dependent

plt.scatter(t, linear_susceptibility.imag)
# limits according to the research paper
plt.xlim(0, 2e-12)
# logarithmic scale according to research paper
plt.yscale('log')
plt.title('X vs Time(s) ')
plt.ylabel('X')
plt.xlabel('Time(s)')
plt.show()


" Separation of Functions "

# Gamma
plt.scatter(q, gamma)
plt.xlim(0, 2e-12)
plt.title('Gamma vs Time(s) ')
plt.ylabel('Gamma')
plt.xlabel('Time(s)')
plt.show()

# n_s
plt.scatter(q, n_s)
plt.xlim(0, 2e-12)
plt.title('n_s vs Time(s) ')
plt.xlabel('Time(s)')
plt.ylabel('n_s')
plt.show()

# np.absolute(np.exp(- 1j * omega_s * t) - 1) ** 2
plt.scatter(t, (np.absolute(np.exp(- 1j * omega_s * t) - 1)) ** 2)
plt.xlim(0, 2e-12)
plt.title('np.absolute(np.exp(- 1j * omega_s * t) - 1) ** 2) (from X_T) vs Time(s) ')
plt.xlabel('Time(s)')
plt.ylabel('np.exp(- 1j * omega_s * t) - 1')
plt.show()

# np.exp(- 1j * omega_s * t) - 1
plt.scatter(t, (np.exp(- 1j * omega_s * t) - 1))
plt.xlim(0, 2e-12)
plt.title('np.exp(- 1j * np.outer(omega_s, t)) - 1) (from X_o) vs Time(s) ')
plt.xlabel('Time(s)')
plt.ylabel('np.exp(- 1j * omega_s * t) - 1')
plt.show()

# np.exp(- 1j * big_omega * t)
plt.plot(t, np.exp(- 1j * big_omega * t))
plt.xlim(0, .25e-12)
plt.title('np.exp(- 1j * big_omega * t (from X) vs Time(s) ')
plt.xlabel('Time(s)')
plt.ylabel('np.exp(- 1j * big_omega * t')
plt.show()


