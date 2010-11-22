#!/usr/bin/env python

from pprint import pprint
import numpy as np
import sys
from numpy import tile
import matplotlib as m
import matplotlib.pyplot as plt
import scipy.io as io

np.set_printoptions(threshold=2000)


def Z(lam, Kw, f, ND, dD):
  """Given Kw, scattering f, ND and dD, computes the reflectivity"""
  return 4*lam**4 / (np.pi**4 * np.abs(Kw)**2) * (np.abs(f)**2 * ND * dD).sum(axis=1) 
  
def get_Zdr(Zh, Zv):
  """Computes Zdr from the individual axis' reflectivity"""
  return 10*np.log10(Zh/Zv)

def get_Zdr_lambda(caplam, D, dD, fa_pi, fb_pi):
  ND =  np.exp(-caplam * D)
  a = (np.abs(fa_pi)**2 * ND * dD).sum(axis=1)
  b = (np.abs(fb_pi)**2 * ND * dD).sum(axis=1)

  return 10*np.log10(a/b)

# get az, zh, abd zdr
file       = "KOUN_20050513-083020.mat"
data       = io.loadmat(file)
data       = data['RawData'][0,0]
Zh         = np.ma.masked_array(data['Zh'], data['Zh'] < 0)
Zdr        = np.ma.masked_array(data['Zdr'], data['Zdr'] < 0)

data = np.load('fwscamp.npz')
fa_pi = data['fapi']
fb_pi = data['fbpi']
fa_0   = data['fa0']
fb_0   = data['fb0']

data = io.loadmat('total.mat')
total = data['total']

# find out where we've got rain
mask = (total != 7) & (total != 8) & (total != 9)
Zh = np.ma.masked_array(Zh,mask)
Zdr = np.ma.masked_array(Zdr,mask)

D,dD = np.linspace(0.08,8,100,retstep=True)
caplam = np.linspace(0,20,100)
caplam = np.tile(caplam, (caplam.shape[0],1))
caplam = np.transpose(caplam)

Zdr_lambda = get_Zdr_lambda(caplam, D, dD, fa_pi, fb_pi)



# Part 2
mu = -0.0201 * caplam**2 + 0.902*caplam - 1.718

R = 1.42e-2 * Zh**0.77 * Zdr**-1.67

#ND = N0*D**mu * np.exp(-caplam * D)


fig = plt.figure(figsize=(15,9));
ax = plt.axes()
#cax = ax.pcolor(R)
ax.plot(caplam[:,0], Zdr_lambda)
#cb = fig.colorbar(cax)
#ax.axis([-150, 150, -150, 150])
ax.set_title("$Z_{dr}$ - $\Lambda$ relation")
ax.set_ylabel(r'$Z_{dr}$ (dB)')
ax.set_xlabel("$\Lambda$ ($mm^{-1}$)")
plt.savefig("zdr-lam.png")
#plt.show()

