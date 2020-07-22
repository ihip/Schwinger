import matplotlib.pyplot as plt

list_x = []
list_y = []
list_dy = []

f = open("pbp.out", "r")
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

plt.title(r"Wilson fermions, $\beta$ = 2.0, L = 16")
plt.xlabel(r"$\kappa$")
plt.ylabel(r"$<\bar{\Psi}\Psi>$")
plt.errorbar(list_x, list_y, yerr = list_dy, linestyle = 'None', marker = '.')
plt.savefig("pbp.svg")

