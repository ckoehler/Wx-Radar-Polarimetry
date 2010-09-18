#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

######## read and massage data #########
file = "hm2p2.txt"
l_D = []
l_sigma_t = []
l_sigma_s = []
l_sigma_b = []


# read the data from the file
for line in open(file):
  cols = line.strip().split("   ")
  l_D.append(cols[0].strip())
  l_sigma_t.append(cols[1].strip())
  l_sigma_s.append(cols[2].strip())
  l_sigma_b.append(cols[3].strip())

# pop off the column labels
l_D.pop(0)
l_sigma_t.pop(0)
l_sigma_s.pop(0)
l_sigma_b.pop(0)

# convert to numpy arrays for the mathy awesome
D = np.array(l_D, dtype=float)
sigma_t = np.array(l_sigma_t, dtype=float)
sigma_s = np.array(l_sigma_s, dtype=float)
sigma_b = np.array(l_sigma_b, dtype=float)
sigma_a = sigma_t - sigma_s
# now we've got all the data we need in a nice set of
# numpy arrays

# these are given
e_r = np.complex(41,41)
e_imag = np.imag(e_r)

# these are easily calculated.
# wavelength in mm, to match D
lam = 30
a = D / 2
k = 2 * np.pi / lam


########## Mie #################
# geometrical cross section:
sigma_g = np.pi * np.power(a,2)

# efficiency factors
Q_t = sigma_t / sigma_g
Q_s = sigma_s / sigma_g
Q_b = sigma_b / sigma_g
Q_a = sigma_a / sigma_g

# scattering albedo
w_0 = sigma_s / sigma_t


########## Rayleigh ##############
# all the cross sections
r_sigma_s = 8 * np.pi * k**4 * a**6 / 3 * np.abs((e_r - 1)/(e_r + 2))**2
r_sigma_b = 4 * np.pi * k**4 * a**6 * np.abs((e_r - 1)/(e_r + 2))**2
r_sigma_a = 4/3 * np.pi * k * a**3 * e_imag * np.abs(3/(e_r + 2))**2
r_sigma_t = r_sigma_a + r_sigma_s

# efficiency factors
r_Q_t = r_sigma_t / sigma_g
r_Q_s = r_sigma_s / sigma_g
r_Q_b = r_sigma_b / sigma_g
r_Q_a = r_sigma_a / sigma_g

# scattering albedo
r_w_0 = r_sigma_s / r_sigma_t

############ Plots #############
# cross sections first
plt.figure(figsize=(15,9));
ax = plt.axes()
ax.plot(D, sigma_t, D, sigma_s, D, sigma_b, D, sigma_a, D, sigma_g, "-.", D, r_sigma_s, "--", D, r_sigma_b, "--", D, r_sigma_a, "--", D, r_sigma_t, "--")
ax.set_title("Cross Sections")
ax.set_xlabel("D (mm)")
ax.set_ylabel("cross sections $\sigma$")
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim((1e-1, 1e3))
ax.set_ylim((1e-3, 1e5))
ax.legend(["Mie $\sigma_t$","Mie $\sigma_s$","Mie $\sigma_b$", "Mie $\sigma_a$", "Geometrical $\sigma_g$", "Rayleigh $\sigma_s$", "Rayleigh $\sigma_b$", "Rayleigh $\sigma_a$", "Rayleigh $\sigma_t$"], loc="lower right")
plt.savefig("crosssections.png")

# then quality factors
plt.figure(figsize=(15,9));
ax = plt.axes()
ax.plot(D, Q_t, D, Q_s, D, Q_b, D, Q_a, D, r_Q_a, '--', D, r_Q_s, '--', D, r_Q_b, '--', D, r_Q_t, '--')
ax.set_xlabel("D (mm)")
ax.set_ylabel("Quality factors")
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim((1e-1, 1e2))
ax.set_ylim((1e-3, 1e1))
ax.set_title("Quality Factors")
ax.legend(["Mie $Q_t$","Mie $Q_s$","Mie $Q_b$", "Mie $Q_a$","Rayleigh $Q_a$","Rayleigh $Q_s$","Rayleigh $Q_b$", "Rayleigh $Q_t$"], loc="upper left")
plt.savefig("qualityfactors.png")

# lastly, albedos
plt.figure(figsize=(15,9));
ax = plt.axes()
ax.plot(D, w_0, D, r_w_0)
ax.set_xlabel("D (mm)")
ax.set_ylabel("Scattering Albedos")
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_title("Scattering Albedo")
ax.legend(["Mie $w_0$","Rayleigh $w_0$"], loc="upper left")
plt.savefig("albedos.png")

plt.show()
