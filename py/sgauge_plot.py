import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt

list_x = []
list_y = []
list_dy = []

f = open("sgauge.out", "r")
for line in f:
    if line[0] == " ":
        data = line.split()
        x = float(data[0])
        y = float(data[1])
        dy = float(data[2])
        print(x, y, dy)
        list_x.append(x)
        list_y.append(y)
        list_dy.append(dy)
f.close()

plt.title(r"pure U(1) HMC, L = 16")
plt.xlabel(r"$\beta$")
plt.ylabel(r"$<S_g>$")

# analytic result
x = np.linspace(0, 8, 100)
plt.ylim([0, 0.7])
plt.plot(x, x * (1.0 - sp.i1(x) / sp.i0(x)), color = 'red', label = 'analytic')
plt.legend()

# numerical data
plt.errorbar(list_x, list_y, yerr = list_dy, linestyle = 'None', marker = '.')

plt.savefig("sgauge-HMC.svg")

