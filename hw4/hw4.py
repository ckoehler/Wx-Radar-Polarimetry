#!/usr/bin/env python

import numpy as np
import scipy.constants as consts
import matplotlib.pyplot as plt


####################
### Functions ######
####################

def axis_ratio(D):
  """calculate axis ratio for raindrops"""
  return 0.995 + 0.0251*D - 0.0364*D**2 + 0.005303*D**3 - 0.0002492*D**4

def gsquared(gamma):
  """find g**2"""
  return (1/gamma**2) - 1

def L_z(D):
  """Get L factor for z-axis"""
  g2 = gsquared(axis_ratio(D))
  g = np.sqrt(g2)
  return (1 + g2)/g2 * (1 - (1/g * np.arctan(g)))

def L_y(lz):
  """computer L factor for y-axis, which depends on L_z in this case"""
  return (1 - lz)/2

def scattering_amplitude(k, D_e, L, e_r):
  """docstring for f_a"""
  a = k**2 * (D_e / 2)**3
  b = e_r - 1
  return a * b / ( 3 * (1 + L * b))


######## read and massage data #########
file = "hm4.dat"
l_D = []
l_fa = []
l_fb = []


# this makes it skip the first line, which is just labels
firstline = True
# read the data from the file
for line in open(file):
  cols = line.strip().split("  ")
  if not firstline:
    l_D.append(cols[0].strip())
    l_fa.append(np.complex(float(cols[1].strip()), float(cols[2].strip())))
    l_fb.append(np.complex(float(cols[3].strip()), float(cols[4].strip())))
  firstline = False
  


# convert to numpy arrays for the mathy awesome
D = np.array(l_D, dtype=float)
t_fa = np.array(l_fa, dtype=complex)
t_fb = np.array(l_fb, dtype=complex)
# now we've got all the data we need in a nice set of
# numpy arrays



############# Given #############
# frequency of the S-band radar
f = 2.87e9
# convert wavelength to mm
lam = consts.c / f * 1000
# and the wave number
k = 2 * np.pi / lam
# effective dielectric constant is also given
e_r = np.complex(80,18)

# generate D from 0.1 to 8mm, 80 data points
D = np.linspace(0.1, 8, 80)

############# Computation ############

# do all the calculations using the methods declared at the top.
gamma = axis_ratio(D)
L_z = L_z(D)
L_y = L_y(L_z)
f_a = scattering_amplitude(k, D, L_y, e_r)
f_b = scattering_amplitude(k, D, L_z, e_r)

############## Plot it all ##############
plt.figure(figsize=(15,9));
ax = plt.axes()
ax.plot(D, np.real(f_a), D, np.real(f_b), D, np.real(t_fa), "--", D, np.real(t_fb), "--")
ax.set_title("Scattering amplitude, real")
ax.set_xlabel("D (mm)")
ax.set_ylabel("Scattering amplitude $\Re(f)$, mm")
ax.set_yscale('log')
ax.legend(["$f_a$", "$f_b$", "T-matrix $f_a$", "T-matrix $f_b$"], loc="upper left" )
plt.savefig("scatteringamps_real.png")

plt.figure(figsize=(15,9));
ax = plt.axes()
ax.plot(D, np.abs(np.imag(f_a)), D, np.abs(np.imag(f_b)), D, np.abs(np.imag(t_fa)), "--", D, np.abs(np.imag(t_fb)), "--")
ax.set_title("Scattering amplitude, imaginary")
ax.set_xlabel("D (mm)")
ax.set_ylabel("Scattering amplitude $|\Im(f)|$, mm")
ax.set_yscale('log')
ax.legend(["$f_a$", "$f_b$", "T-matrix $f_a$", "T-matrix $f_b$"], loc="upper left" )
plt.savefig("scatteringamps_imag.png")

plt.show()
