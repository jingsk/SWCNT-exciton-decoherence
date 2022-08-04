import numpy as np
import warnings
import matplotlib.pyplot as plt
from scipy.integrate import quad

warnings.filterwarnings('ignore')

" CONSTANTS "

# Boltzmann constant
k_b = 1.380649e-23  # J/k
# h bar
h_bar = 1.054571817e-34  # Js
# temperature
T = 6  # K
# linear mass density
p = 1.67e-15  # kg/m
# vs
v = 19900  # m/s
# deformation potential
Ds = 2.243047191e-18  # J
# length of nanotube
L = 1.883e-9  # m
# confinement length
sigma = 3e-9  # m
# energy of the excited state
omega_constant = 1.41 * 1.60218e-19 / h_bar  # Hz

# wave vector
q_array = (np.linspace(0, 0.2, num=1000) + .01) * 1 / L  # 1/m

# time
t_array = np.linspace(0, 3, num=1000) * 1e-12  # s


# dimensionless coupling strength
def gamma(q):
    # deformation potential couplings
    G = Ds * q / np.sqrt(2 * p * L * h_bar * v * q)
    # form factor
    F = np.exp(-(q ** 2 * sigma ** 2) / 4)
    # linear dispersion
    w = v * q
    return G * F / w


# phonon occupation number
def n(q):
    return (np.exp((h_bar * v * q) / (k_b * T)) - 1) ** -1


# linear susceptibility temperature independent
def X_o(q, t):
    return gamma(q) ** 2 * (1 - np.cos(v * q * t))


exponents = - (L / np.pi) * np.array([quad(X_o, 0, np.pi / L, args=t)[0] for t in t_array])
norm_X_o = np.exp(exponents)


# linear susceptibility temperature dependent
def X_T(q, t):
    return n(q) * gamma(q) ** 2 * (2 - 2 * np.cos(v * q * t))


exponents = - (L / np.pi) * np.array([quad(X_T, 0, np.pi / L, args=t)[0] for t in t_array])
norm_X_T = np.exp(exponents)


plt.plot(t_array, norm_X_o, t_array, norm_X_T)
plt.title('|X_T| and |X_o| vs Time (ps)')
plt.yscale('log')
plt.legend(['X_o', 'X_T'])
plt.ylim([0.08, 1.1])
plt.ylabel('|X_T| and |X_o|')
plt.xlabel('Time(ps)')
plt.show()


# linear susceptibility
def X(q, t):
    return gamma(q) ** 2 * (1 + 2 * n(q)) * (1 - np.cos(v * q * t))


# high temperature
T = 6  # K
exponents = - (L / np.pi) * np.array([quad(X, 0, np.pi / L, args=t)[0] for t in t_array])
norm1 = np.exp(exponents)

# medium temperature
T = 3  # K
exponents = - (L / np.pi) * np.array([quad(X, 0, np.pi / L, args=t)[0] for t in t_array])
norm2 = np.exp(exponents)

# low temperature
T = .01  # K
exponents = - (L / np.pi) * np.array([quad(X, 0, np.pi / L, args=t)[0] for t in t_array])
norm4 = np.exp(exponents)

# exponential slope on logarithmic scale
y = 1.4 * t_array ** .955

# plotting
plt.plot(t_array, norm1, t_array, norm2, t_array, norm4, t_array, y)
plt.title('Linear Susceptibility (Research Paper) vs Time (ps)')
plt.yscale('log')
plt.legend(['6 K', '3 K', '0.01 K'])
plt.ylim([0.08, 1.1])
plt.ylabel('|X(t)|')
plt.xlabel('Time(ps)')
plt.show()

# (5,6) SWNTs chirality

sigma = 3.3e-10  # m
T = 6  # K

exponents = - (L / np.pi) * np.array([quad(X, 0, np.pi / L, args=t)[0] for t in t_array])
norm5 = np.exp(exponents)

# medium temperature
T = 3  # K
exponents = - (L / np.pi) * np.array([quad(X, 0, np.pi / L, args=t)[0] for t in t_array])
norm6 = np.exp(exponents)

# low temperature
T = .01  # K
exponents = - (L / np.pi) * np.array([quad(X, 0, np.pi / L, args=t)[0] for t in t_array])
norm7 = np.exp(exponents)

# plotting
plt.plot(t_array, norm5, t_array, norm6, t_array, norm7, t_array, y)
plt.title('Linear Susceptibility ((5,6) SWNTs) vs Time (ps)')
plt.yscale('log')
plt.legend(['6 K', '3 K', '0.01 K'])
plt.ylim([0.08, 1.1])
plt.ylabel('|X(t)|')
plt.xlabel('Time(ps)')
plt.show()



