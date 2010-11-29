#!/usr/bin/env python

from pprint import pprint
import numpy as np
import sys
from numpy import tile
import matplotlib as m
import scipy.constants as consts
import matplotlib.pyplot as plt
import scipy.io as io
import scipy.interpolate as ip

np.set_printoptions(threshold=2000)

def plot(x,y,title,xlabel,ylabel,figname):
  """plot stuff"""
  fig = plt.figure(figsize=(15,9));
  ax = plt.axes()
  ax.plot(x, y)
  ax.set_title(title)
  ax.set_xlabel(xlabel)
  ax.set_ylabel(ylabel)
  plt.savefig(figname)

def hist(x,n,title,xlabel,ylabel,figname, axis=[]):
  """plot stuff"""
  fig = plt.figure(figsize=(15,9));
  ax = plt.axes()
  ax.hist(x, n)
  ax.set_title(title)
  ax.set_xlabel(xlabel)
  if(axis):
    ax.axis(axis)
  ax.set_ylabel(ylabel)
  plt.savefig(figname)

def Z(lam, Kw, f, ND, dD):
  """Given Kw, scattering f, ND and dD, computes the reflectivity"""
  return 4*lam**4 / (np.pi**4 * np.abs(Kw)**2) * (np.abs(f)**2 * ND * dD).sum(axis=1) 
  
def get_Zdr(Zh, Zv):
  """Computes Zdr from the individual axis' reflectivity"""
  return 10*np.log10(Zh/Zv)

def get_Zdr_lambda(ND, dD, fa_pi, fb_pi):
  a = (np.abs(fa_pi)**2 * ND * dD).sum(axis=1)
  b = (np.abs(fb_pi)**2 * ND * dD).sum(axis=1)

  return 10*np.log10(a/b)

def get_N0(Zh,Kw,lam,f,ND,dD):
  """compute N0 from Zh, given Lambda"""
  a = np.pi**4 * np.abs(Kw)**2 / (4*lam**4)
  b = (np.abs(f)**2 * ND * dD).sum(axis=0)
  c = Zh * a
  d = b**-1
  return c * d

def get_N0_cg(Zh,Kw,lam,f,ND,dD):
  """compute N0 from Zh, given Lambda"""
  a = np.pi**4 * np.abs(Kw)**2 / (4*lam**4)
  b = (np.abs(f)**2 * ND * dD).sum(axis=0)
  c = Zh * a
  d = b**-1
  return c * d

def nthMoment(D,Nd, dD,n):
  """docstring for nthMoment"""
  return ((D**n) * Nd * dD).sum(axis=0)

def RainRate(D, Nd, Vd, dD) :
    # units are mm*hr^-1
    return 6e-4 * np.pi * nthMoment(D, Vd*Nd, dD, 3)

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

f = 2.8e9
c = consts.c
# convert lambda to mm
lam = c/f * 1e3

er = 80.2093369315+17.1572908221j
Kw = (er - 1)/(er + 2)

# find out where we've got rain
mask     = (total != 7) & (total != 8) & (total != 9)
Zh       = np.ma.masked_array(Zh,mask)
Zdr      = np.ma.masked_array(Zdr,mask)

D,dD     = np.linspace(0.08,8,100,retstep=True)
caplam   = np.linspace(0,20,100)
caplam   = np.tile(caplam, (caplam.shape[0],1))
caplam   = np.transpose(caplam)
caplam1d = caplam[:,0]

# exp model first
ND =  np.exp(-caplam * D)
Zdr_lambda = get_Zdr_lambda(ND, dD, fa_pi, fb_pi)
sort_order = np.argsort(Zdr_lambda)

fit = ip.InterpolatedUnivariateSpline(Zdr_lambda[sort_order], caplam1d[sort_order])
fitted_Zdr = fit(np.abs(Zdr.flatten()))[::-1]
all_caplam = np.ma.masked_array(fitted_Zdr, fitted_Zdr <= -10) 

ND = np.exp(-caplam1d * D)
N0 = get_N0(Zh.flatten(), Kw, lam, fa_pi, ND, dD)
N0 = np.ma.masked_array(N0, N0 <= -100)

#plot(caplam1d, Zdr_lambda,"$Z_{dr}$ - $\Lambda$ relation, Exp model",ylabel='$Z_{dr}$ (dB)',xlabel="$\Lambda$ ($mm^{-1}$)",figname="zdr-lam-exp.png")
#hist(all_caplam.compressed(),1000, "$\Lambda$ for all $Z_{dr}$", "$\Lambda$", "Number", "all-caplam.png", axis=[-5,10,0,6500])
#hist(N0.compressed(),300, "$N_0$ for all $Z_{dr}$", "$N_0$", "Number", "N0.png")

# now C-G ND
mu = -0.0201 * caplam1d**2 + 0.902*caplam1d - 1.718
ND = D**mu * np.exp(-caplam * D)
Zdr_lambda = get_Zdr_lambda(ND, dD, fa_pi, fb_pi)
sort_order = np.argsort(Zdr_lambda)

fit = ip.InterpolatedUnivariateSpline(Zdr_lambda[sort_order], caplam1d[sort_order])
fitted_Zdr = fit(np.abs(Zdr.flatten()))[::-1]
all_caplam = np.ma.masked_array(fitted_Zdr, fitted_Zdr <= -10) 

ND = D**mu * np.exp(-caplam1d * D)
N0 = get_N0_cg(Zh.flatten(), Kw, lam, fa_pi, ND, dD)
N0 = np.ma.masked_array(N0, N0 <= -100)

#plot(caplam1d, Zdr_lambda, "$Z_{dr}$ - $\Lambda$ relation, C-G model",ylabel='$Z_{dr}$ (dB)', xlabel="$\Lambda$ ($mm^{-1}$)",figname="zdr-lam-cg.png")
#hist(all_caplam.compressed(),1000, "$\Lambda$ for all $Z_{dr}$", "$\Lambda$", "Number", "all-caplam-cg.png", axis=[2,8,0,3000])
#hist(N0.compressed(),300, "$N_0$ for all $Z_{dr}$", "$N_0$", "Number", "N0-cg.png")


# Rainfall rate and D_m
v = -0.1021 + 4.932*D- 0.9551*D**2 + 0.07934*D**3 - 0.002362*D**4
D2d   = np.tile(D, (all_caplam.shape[0],1))
#mu   = np.tile(mu, (all_caplam.shape[0],1))

v   = np.transpose(np.tile(v, (all_caplam.shape[0],1)))

#ND_ret = N0 * np.transpose(D2d**mu) * np.exp(-all_caplam * np.transpose(D))
#print ND_ret.shape
#R_calc = RainRate(np.transpose(D2d), ND_ret, v, dD)

#print R_calc.shape
#ND = N0*D2d**mu * np.exp(-caplam * D2d)


### PART 2 ######

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


fitted_Zdr = fit(np.abs(koun_Zdr.flatten()))[::-1]
all_caplam = np.ma.masked_array(fitted_Zdr, fitted_Zdr <= -10) 

ND = D**mu * np.exp(-caplam1d * D)
N0 = get_N0_cg(koun_Zh.flatten(), Kw, lam, fa_pi, ND, dD)
N0 = np.ma.masked_array(N0, N0 <= -100)

hist(all_caplam.compressed(),1000, "$\Lambda$ for all $Z_{dr}$", "$\Lambda$", "Number", "2-Lambda.png")#, axis=[2,8,0,3000])
hist(N0.compressed(),300, "$N_0$ for all $Z_{dr}$", "$N_0$", "Number", "2-N0.png")

# DSD data
# get measured DSD data
file = "dsddata_20050513.mat"
data = io.loadmat(file)

# define old and new length x values to
# interpolate over next.
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

# calculate v and R with given data
dD_meas = 0.2
v = -0.1021 + 4.932*D_meas - 0.9551*D_meas**2 + 0.07934*D_meas**3 - 0.002362*D_meas**4
R_calc = 6e-4 * np.pi * (D_meas**3 * v * ND_meas * dD_meas).sum(axis=1)

ND = D**mu * np.exp(-caplam * D)
D = np.tile(D, (D.shape[0],1))
print D_meas.shape
print ND_meas.shape
print D.shape
print ND.shape
Dt_meas = nthMoment(D_meas, ND_meas, dD_meas, 0)
Dt_retr = nthMoment(D, ND, dD, 0)
print Dt_meas.shape
print Dt_retr.shape


fig = plt.figure(figsize=(15,9));
ax = plt.axes()
ax.plot(xx, Dt_meas)
ax.plot(xx, Dt_retr)
ax.legend(["meas", 'retr'])
ax.set_title("Dt")
ax.set_xlabel("t")
ax.set_ylabel("Dt")
plt.savefig("2-Dt.png")