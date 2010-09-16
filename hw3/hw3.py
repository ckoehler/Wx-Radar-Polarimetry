#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

# read and massage data
file = "hm2p2.txt"
l_D = []
l_sigma_t = []
l_sigma_s = []
l_sigma_b = []


for line in open(file):
  cols = line.strip().split("   ")
  l_D.append(cols[0].strip())
  l_sigma_t.append(cols[1].strip())
  l_sigma_s.append(cols[2].strip())
  l_sigma_b.append(cols[3].strip())

l_D.pop(0)
l_sigma_t.pop(0)
l_sigma_s.pop(0)
l_sigma_b.pop(0)

D = np.array(l_D, dtype=float)
sigma_t = np.array(l_sigma_t, dtype=float)
sigma_s = np.array(l_sigma_s, dtype=float)
sigma_b = np.array(l_sigma_b, dtype=float)
# now we've got all the data we need in a nice
# numpy array

# these are given
e = np.complex(41,41)
lam = 0.03


# geometrical cross section:
sigma_g = np.pi * np.power(D,2)

# efficiency factors
Q_t = sigma_t / sigma_g;
Q_s = sigma_s / sigma_g;
Q_b = sigma_b / sigma_g;


# plot it all
fig = plt.figure()
ax = fig.add_subplot(2,1,1)
ax.plot(D, sigma_t, D, sigma_s, D, sigma_b, D, sigma_g)
ax.set_title("Cross Sections")
ax.set_xlabel("D (mm)")
ax.legend(["Mie $\sigma_t$","Mie $\sigma_s$","Mie $\sigma_b$", "Geometrical $\sigma_g$"], loc="upper left")

ax = fig.add_subplot(2,1,2)
ax.plot(D, Q_t, D, Q_s, D, Q_b)
ax.set_xlabel("D (mm)")
ax.set_title("Quality Factors")
ax.legend(["Mie $Q_t$","Mie $Q_s$","Mie $Q_b$"], loc="upper right")

plt.show()
