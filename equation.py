import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Boltzmann constant (J/K)
k = 1.380649 ** 0.00000000000000000000010
# h bar (J-s)
h = 1.054571817 ** .0000000000000000000000000000000010
# temperature (K)
T = 4
# e
e = 2.718281828459045

# exciton coupling matrix elements (default values)
g = np.random.randint(1, 10, size=(4, 6))


# linear dispersion (default values)
w = np.random.randint(1, 10, size=(1, 6))

# dimensionless coupling strength
gamma = np.divide(g, w)

# phonon occupation number
n = []

for i in w[0]:
    value = (e ** ((h * i)/(k * T)) - 1) ** -1
    n.append(value)


print(g)
print(w)
print(gamma)
print(n)

# Here a sample on how to read a file
# sample_data = pd.read_csv('sample_data.csv')

