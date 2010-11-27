#!/usr/bin/env python

from pprint import pprint
import numpy as np
import sys
from numpy import tile
import matplotlib as m
import matplotlib.pyplot as plt
import scipy.io as io
import scipy.interpolate as ip

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

def get_N0(Zh,Kw,lam,f,caplam,D,dD):
  """compute N0 from Zh, given Lambda"""
  a = np.pi**4 * np.abs(Kw)**2 / (4*lam**4)
  b = (np.abs(f)**2 * np.exp(-caplam * D) * dD).sum(axis=1)
  return Zh * a * b**-1 

def nthMoment(D,ND, dD,n):
  """docstring for nthMoment"""
  return ((D**n) * Nd * dD).sum(axis=0)

def RainRate(D, Nd, Vd) :
    # units are mm*hr^-1
    return 6e-4 * np.pi * nthMoment(D, Vd*Nd, 3)

def Dm(D,ND,dD):
  """volume weighted diamter"""
  return nthMoment(D,ND,dD,4) / nthMoment(D,ND,dD,3)

# get az, zh, abd zdr
file       = "KOUN_20050513-083020.mat"
data       = io.loadmat(file)
data       = data['RawData'][0,0]
Zh_cond    = (data['Zh'] == -99900.0) | (data['Zh'] == -32768.0)
Zdr_cond   = (data['Zdr'] == -99900) | (data['Zdr'] == -32768) | (data['Zdr'] > 4.7) | (data['Zdr'] < -4.7)
Zh         = np.ma.masked_array(data['Zh'], Zh_cond )
Zdr        = np.ma.masked_array(data['Zdr'], Zdr_cond)
data  = np.load('fwscamp.npz')
fa_pi = data['fapi']
fb_pi = data['fbpi']
fa_0  = data['fa0']
fb_0  = data['fb0']

data  = io.loadmat('total.mat')
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

fit = ip.interp1d(Zdr_lambda, caplam[:,0])
#all_caplam = fit(Zdr.compressed())
#print fit(1)
#print fit(4)
print Zdr_lambda.min()
print Zdr_lambda.max()

# Part 2
mu = -0.0201 * caplam**2 + 0.902*caplam - 1.718

R = 1.42e-2 * Zh**0.77 * Zdr**-1.67

#ND = N0*D**mu * np.exp(-caplam * D)

fig = plt.figure(figsize=(15,9));
ax = plt.axes()
ax.plot(Zdr_lambda, caplam[:,0])
ax.plot(np.linspace(0,5,100), fit(np.linspace(0,5,100)))
ax.set_title("$Z_{dr}$ - $\Lambda$ relation")
ax.set_xlabel(r'$Z_{dr}$ (dB)')
ax.set_ylabel("$\Lambda$ ($mm^{-1}$)")
plt.savefig("zdr-lam.png")

fig = plt.figure(figsize=(15,9));
ax = plt.axes()
ax.hist(all_caplam, 100)
#cb = fig.colorbar(cax)
ax.set_title("$\Lambda$ for all $Z_{dr}$")
#ax.set_xlabel(r'$Z_{dr}$ (dB)')
#ax.set_ylabel("$\Lambda$ ($mm^{-1}$)")
plt.savefig("all-caplam.png")


#plt.show()

