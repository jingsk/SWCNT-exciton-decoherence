import numpy as np
import warnings
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')


" CONSTANTS "


# move constant to a different files

# Boltzmann constant (J/K)
K_b = 1.380649e-23
# h bar (J-s)
h_bar = 1.054571817e-34
# temperature (K)
T = 6
# rho (kg/m)
p = 1.67e-15
# vs (m/s)
v = 19900
# deformation potential (J)
Ds = 1.9228e-18
# length of nanotube (m), change
L = 18.83e-10
# confinement length (m), change
sigma = 1e-8
# energy of the excited state (Hz)
omega_constant = 1.41 * 1.60218e-19 / h_bar
# Unknown constant from directly proportional relation
C = 1


" FUNCTIONS "

# wave vector
q = (np.linspace(0, 0.2, num=10) + .01) * 1/L

# time (s)
t = np.linspace(0, 3, num=20000) * 1e-10

# linear dispersion
omega_s = v * q

# deformation potential couplings
G = Ds * q / np.sqrt(2 * p * L * h_bar * omega_s)

# form factor
F = np.exp(-(q ** 2 * sigma ** 2)/4)

# exciton coupling matrix elements
g = G * F

# dimensionless coupling strength
gamma = g / omega_s

# phonon occupation number
n_s = (np.exp((h_bar * omega_s)/(K_b * T)) - 1) ** -1

# temperature_dependent
temperature_dependent = 1j * np.exp(np.absolute(gamma) ** 2 * - n_s
                                    @ (np.absolute(np.exp(- 1j * np.outer(omega_s, t)) - 1)) ** 2)

# temperature_independent
temperature_independent = 1j * np.exp(np.absolute(gamma) ** 2
                                      @ (np.exp(- 1j * np.outer(omega_s, t)) - 1))

# phonon shifted transition frequency
big_omega = omega_constant * h_bar - np.absolute(gamma) ** 2 @ omega_s

# linear susceptibility
linear_susceptibility = C * - 1j * np.exp(- 1j * big_omega * t) \
                        * temperature_independent * temperature_dependent

# norm_Xo_XT
norm_Xo_XT = np.exp( - (np.absolute(gamma) ** 2 * (1 + 2 * n_s) @ (1 - np.cos(np.outer(omega_s, t)))))

# norm_Xo
norm_Xo = np.exp( - (np.absolute(gamma) ** 2 @ (1 - np.cos(np.outer(omega_s, t)))))

# norm_XT
norm_XT = np.exp( - (np.absolute(gamma) ** 2 * n_s @ (2 - 2 * np.cos(np.outer(omega_s, t)))))

plt.scatter(t, norm_Xo_XT, label = "norm_Xo_XT (Analytical)")
plt.scatter(t, norm_Xo, label = "norm_Xo (Analytical)")
plt.scatter(t, norm_XT, label = "norm_XT (Analytical)")

plt.xlim(0, .5e-12)
plt.ylim(0.08, 1)
plt.yscale('log')
plt.title('Norms vs Time(ps) ')
plt.ylabel('Norms')
plt.xlabel('Time(ps)')
plt.legend()
plt.show()

plt.scatter(t, np.absolute(linear_susceptibility), label = "norm_Xo_XT (Coding)")
plt.scatter(t, np.absolute(temperature_independent), label = "norm_Xo (Coding)")
plt.scatter(t, np.absolute(temperature_dependent), label = "norm_XT (Coding)")

plt.xlim(0, .5e-12)
plt.ylim(0.08, 1)
plt.yscale('log')
plt.title('Norms vs Time(ps) ')
plt.ylabel('Norms')
plt.xlabel('Time(ps)')
plt.legend()
plt.show()
