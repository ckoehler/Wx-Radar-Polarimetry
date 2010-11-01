#!/usr/bin/env python

import numpy as np
import scipy.constants as consts
import matplotlib.pyplot as plt
import scipy.io as io



def Z(lam, Kw, f, ND, dD):
  """Given Kw, scattering f, ND and dD, computes the reflectivity"""
  return 4*lam**4 / (np.pi**4 * np.abs(Kw)**2) * (np.abs(f)**2 * ND * dD).sum(axis=1)
  
def Zdr(Zh, Zv):
  """Computes Zdr from the individual axis' reflectivity"""
  return 10*np.log(Zh/Zv)

Ah = 8.686*lam * (np.imag(fa_0) * ND * dD).sum(axis=0)
Av = 8.686*lam * (np.imag(fb_0) * ND * dD).sum(axis=0)

Adp = Ah - Av

Kdp = 180*lam / np.pi * (np.real(dD*fa_0 - dD*fb_0) * ND * dD).sum(axis=1)


######## read and massage data #########
file = "koun.dat"
l_t = []
l_Zh = []
l_Zdr = []


# read the data from the file
for line in open(file):
  cols = line.strip().split("   ")
  l_t.append(cols[0].strip())
  l_Zh.append(cols[1].strip())
  l_Zdr.append(cols[2].strip())

# pop off the column labels
l_t.pop(0)
l_Zh.pop(0)
l_Zdr.pop(0)

# convert to numpy arrays for the mathy awesome
koun_t = np.array(l_t, dtype=float)
koun_Zh = np.array(l_Zh, dtype=float)
koun_Zdr = np.array(l_Zdr, dtype=float)


# repeat for second file
file = "SCTT_RAIN_fw100.dat"
l_size = []
l_fa_180 = []
l_fb_180 = []
l_fa_0 = []
l_fb_0 = []

# this makes it skip the first line, which is just labels
firstline = True
# read the data from the file
for line in open(file):
  cols = line.strip().split("  ")
  if not firstline:
    l_size.append(cols[0].strip())
    l_fa_180.append(np.complex(float(cols[1].strip()), float(cols[2].strip())))
    l_fb_180.append(np.complex(float(cols[3].strip()), float(cols[4].strip())))
    l_fa_0.append(np.complex(float(cols[5].strip()), float(cols[6].strip())))
    l_fb_0.append(np.complex(float(cols[7].strip()), float(cols[8].strip())))
  firstline = False

# convert to numpy arrays for the mathy awesome
size = np.array(l_size, dtype=float)
fa_180 = np.array(l_fa_180, dtype=complex)
fb_180 = np.array(l_fb_180, dtype=complex)
fa_0 = np.array(l_fa_0, dtype=complex)
fb_0 = np.array(l_fb_0, dtype=complex)

# now we've got all the data we need in a nice set of
# numpy arrays





###### givens #######
f = 2.8e9
c = consts.c

# convert lambda to mm
lam = c/f * 1e3
#t2 = np.linspace(np.min(t),np.max(t),100)


##### PART 1 #######

##### calculations #####
R = np.linspace(1,100,100)
R = np.tile(R, (R.shape[0],1))

# R is in mm/hr
R = np.transpose(R)
delta = 4.1*R**(-0.21)

# this D is identical to `size` above
# D is in mm
D,dD = np.linspace(0.08,8,100,retstep=True)

# each row corresponds to one set of ND as a function of one R
ND = 8000*np.exp(-delta*D)

# dielectric constant of water at 10 degrees C
Kw = 80.2093369315+17.1572908221j
Zh = 4*lam**4 / (np.pi**4 * np.abs(Kw)**2) * (np.abs(fa_0)**2 * ND * dD).sum(axis=1)
Zh = Z(lam, Kw, fa_0, ND, dD)
Zv = Z(lam, Kw, fb_0, ND, dD)
Zv = 4*lam**4 / (np.pi**4 * np.abs(Kw)**2) * (np.abs(fb_0)**2 * ND * dD).sum(axis=1)
Zdr = 10*np.log(Zh/Zv)

Ah = 8.686*lam * (np.imag(fa_0) * ND * dD).sum(axis=0)
Av = 8.686*lam * (np.imag(fb_0) * ND * dD).sum(axis=0)

Adp = Ah - Av

Kdp = 180*lam / np.pi * (np.real(dD*fa_0 - dD*fb_0) * ND * dD).sum(axis=1)

#print Ah.shape
#print Av.shape
#print Adp.shape
#print Kdp.shape

fig = plt.figure(figsize=(15,9));
ax = fig.add_subplot(2,1,1)
ax.plot(R[:,0], Zh, R[:, 0], Zv)
ax.set_title("$Z_h$ and $Z_v$")
ax.set_ylabel(r'$mm^6/m^3$')
ax.legend(["$Z_h$","$Z_v$"], loc="upper left")

ax = fig.add_subplot(2,1,2)
ax.plot(R[:,0], Zdr)
ax.set_title("$Z_{dr}$")
ax.set_xlabel("rainfall rate R (mm/hr)")
ax.set_ylabel(r'dB')
plt.savefig("1-zhzv.png")


fig = plt.figure(figsize=(15,9));
ax = fig.add_subplot(2,1,1)
ax.plot(R[:,0], Ah, R[:, 0], Av)
ax.set_title("$A_h$ and $A_v$")
ax.set_ylabel(r'$dB/km$')
ax.legend(["$A_h$","$A_v$"], loc="upper left")

ax = fig.add_subplot(2,1,2)
ax.plot(R[:,0], Adp)
ax.set_title("$A_{dp}$")
ax.set_xlabel("rainfall rate R (mm/hr)")
ax.set_ylabel(r'dB')
plt.savefig("1-AhAv.png")


fig = plt.figure(figsize=(15,9));
ax = plt.axes()
ax.plot(R[:,0], Kdp)
ax.set_title("$K_{dp}$")
ax.set_ylabel(r'$deg/km$')
ax.set_xlabel("rainfall rate R (mm/hr)")
plt.savefig("1-Kdp.png")



##### PART 2 #######

# get measured DSD data
file = "dsddata_20050513.mat"
data = io.loadmat(file)

D_meas = data['dsd_data'][:, 2, :]
Nd_meas = data['dsd_data'][:, 5, :]
dD_meas = 0.2


v = -0.1021 + 4.932*D_meas - 0.9551*D_meas**2 + 0.07934*D_meas**3 - 0.002362*D_meas**4

R_calc = 6e-4 * np.pi * (D_meas**3 * v * ND_meas * dD_meas).sum(axis=0)


##### PART 3 #######
#plt.figure(figsize=(15,9));
#ax = plt.add_subplot(1,4,1)
##ax.plot(R, Zh, R, Zv, R, Ah, R, Av, R, Adp, R, Kdp )
#ax.plot(R, Zh, R, Zv)
##ax.plot(D, ND[45])
#ax.set_title("Comparison of measured vs calculated $Z_hh$")
#ax.set_xlabel("rainfall rate R")
##ax.set_yscale('log')
#ax.set_ylabel("Zh")
#ax.legend(["$Z_h$","$Z_v$", "$A_h$", "$A_v$", "$A_dp$", "$K_dp$"], loc="upper right")
#plt.savefig("comparison.png")
#plt.show()
