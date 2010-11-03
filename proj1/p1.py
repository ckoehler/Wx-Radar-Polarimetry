#!/usr/bin/env python

import numpy as np
import scipy.constants as consts
import matplotlib.pyplot as plt
import scipy.io as io
import scipy.interpolate as ip



def Z(lam, Kw, f, ND, dD):
  """Given Kw, scattering f, ND and dD, computes the reflectivity"""
  return 4*lam**4 / (np.pi**4 * np.abs(Kw)**2) * (np.abs(f)**2 * ND * dD).sum(axis=1)
  
def get_Zdr(Zh, Zv):
  """Computes Zdr from the individual axis' reflectivity"""
  return 10*np.log(Zh/Zv)

def A(lam, f, ND, dD):
  """Get attenuation for given parameters"""
  print "ND"
  print ND.shape
  return 8.686*lam * (np.imag(f) * ND * dD).sum(axis=1) * 1e-3

def get_Adp(Ah, Av):
  """Differential attenuation"""
  return Ah - Av

def get_Kdp(lam, fa, fb, ND, dD):
  """Differential phase"""
  return 180*lam / np.pi * (np.real(fa - fb) * ND * dD).sum(axis=1) * 1e-3


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

# dielectric constant of water at 10 degrees C
Kw = 80.2093369315+17.1572908221j

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

Zh = Z(lam, Kw, fa_180, ND, dD)
Zv = Z(lam, Kw, fb_180, ND, dD)
Zdr = get_Zdr(Zh, Zv)
Ah = A(lam, fa_0, ND, dD)
Av = A(lam, fb_0, ND, dD)
Adp = get_Adp(Ah, Av)
Kdp = get_Kdp(lam, fa_0, fb_0, ND, dD)

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
ax.legend(["$A_h$","$A_v$"], loc="upper right")

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

# this is time data, so make a new time array
x = np.linspace(1,41,41)
xx = np.linspace(1,41,100)

# interpolate the data to fit into 100 data points
D_meas = np.transpose(data['dsd_data'][:, 2, :])
s = ip.interp1d(x, D_meas, axis=1, kind='cubic')
D_meas = s(xx)

# interpolate the data to fit into 100 data points
ND_meas = np.transpose(data['dsd_data'][:, 5, :])
s = ip.interp1d(x, ND_meas, axis=1, kind='cubic')
ND_meas = s(xx)

dD_meas = 0.2
v = -0.1021 + 4.932*D_meas - 0.9551*D_meas**2 + 0.07934*D_meas**3 - 0.002362*D_meas**4
R_calc = 6e-4 * np.pi * (D_meas**3 * v * ND_meas * dD_meas).sum(axis=1)

meas_Zh = Z(lam, Kw, fa_180, ND_meas, dD_meas)
meas_Zv = Z(lam, Kw, fb_180, ND_meas, dD_meas)
meas_Zv = np.ma.masked_array(meas_Zv, meas_Zv == 0.0)
meas_Zdr = get_Zdr(meas_Zh, meas_Zv)
meas_Ah = A(lam, fa_0, ND_meas, dD_meas)
meas_Av = A(lam, fb_0, ND_meas, dD_meas)
meas_Adp = get_Adp(meas_Ah, meas_Av)
meas_Kdp = get_Kdp(lam, fa_0, fb_0, ND_meas, dD_meas)

t = np.linspace(0,100,481)
fig = plt.figure(figsize=(15,9));
ax = fig.add_subplot(2,1,1)
ax.plot(t, meas_Zh, t, meas_Zv)
ax.set_title("$Z_h$ and $Z_v$")
ax.set_ylabel(r'$mm^6/m^3$')
ax.legend(["$Z_h$","$Z_v$"], loc="upper right")

ax = fig.add_subplot(2,1,2)
ax.plot(t, meas_Zdr)
ax.set_title("$Z_{dr}$")
ax.set_xlabel("Time")
ax.set_ylabel(r'dB')
plt.savefig("2-zhzv.png")


fig = plt.figure(figsize=(15,9));
ax = fig.add_subplot(2,1,1)
ax.plot(t, meas_Ah, t, meas_Av)
ax.set_title("$A_h$ and $A_v$")
ax.set_ylabel(r'$dB/km$')
ax.legend(["$A_h$","$A_v$"], loc="upper left")

ax = fig.add_subplot(2,1,2)
ax.plot(t, meas_Adp)
ax.set_title("$A_{dp}$")
ax.set_xlabel("Time")
ax.set_ylabel(r'dB')
plt.savefig("2-AhAv.png")


fig = plt.figure(figsize=(15,9));
ax = plt.axes()
ax.plot(meas_Kdp)
ax.set_title("$K_{dp}$")
ax.set_ylabel(r'$deg/km$')
ax.set_xlabel("Time")
plt.savefig("2-Kdp.png")

fig = plt.figure(figsize=(15,9));
ax = plt.axes()
ax.plot(R_calc)
ax.set_title("Rainfall rate $R$")
ax.set_ylabel(r'$mm/hr$')
ax.set_xlabel("Time")
plt.savefig("2-R.png")



## associate rainfall rate and reflectivity
fig = plt.figure(figsize=(15,9));
ax = fig.add_subplot(2,1,1)
ax.plot(Zh)
ax.scatter(R_calc, meas_Zh)
ax.set_title("$Z_h$ as a function of R")
ax.set_ylabel(r'$Z_h$')
ax.set_xlabel(r'$R$')
ax.legend([r'$Z_h$ from MP', r'$Z_h$ from DSD data'], loc='lower right')
#ax.axis([0,160,0,80])
ax = fig.add_subplot(2,1,2)
ax.plot(Zdr)
ax.scatter(R_calc, meas_Zdr)
ax.set_title("$Z_{dr}$ as a function of R")
ax.set_ylabel(r'$Z_{dr}$')
ax.set_xlabel(r'$R$')
ax.legend([r'$Z_{dr}$ from MP', r'$Z_{dr}$ from DSD data'], loc='lower right')
#ax.axis([0,160,0,10])
plt.savefig("2-zhzdrofr.png")


##### PART 3 #######

# now we just have minutes since 7:00UTC here, same as the DSD data
koun_t = np.array(np.round((koun_t - 7) * 60), dtype='int');

comp_meas_Zh = meas_Zh[koun_t]
comp_meas_Zdr = meas_Zdr[koun_t]



fig = plt.figure(figsize=(15,9));
ax = fig.add_subplot(2,1,1)
ax.plot(koun_Zh)
ax.plot(comp_meas_Zh)
ax.set_title("$Z_h$ from KOUN")
ax.set_xlabel("Time")
ax.set_ylabel(r'dBZ')
ax.legend(['KOUN', 'DSD'], loc='upper right')

ax = fig.add_subplot(2,1,2)
ax.plot(koun_Zdr)
ax.plot(comp_meas_Zdr)
ax.set_title("$Z_{dr}$ from KOUN")
ax.set_xlabel("Time")
ax.set_ylabel(r'dB')
ax.legend(['KOUN', 'DSD'], loc='upper right')
plt.savefig("3-koun.png")
#plt.show()
