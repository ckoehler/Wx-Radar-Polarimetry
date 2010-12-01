#!/usr/bin/env python

from pprint import pprint
import numpy as np
import sys
import os.path
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
  b = (np.abs(f)**2 * ND * dD).sum(axis=1)
  c = Zh * a
  d = b**-1
  return c * d

def get_N02(Zh,Kw,lam,f,ND,dD):
  """compute N0 from Zh, given Lambda"""
  a = np.pi**4 * np.abs(Kw)**2 / (4*lam**4)
  b = (np.abs(f)**2 * ND * dD).sum(axis=0)
  c = Zh * a
  d = 1/b
  return c * d

def nthMoment(D,Nd, dD,n):
  """calculate the nth-moment of the DSD"""
  return ((D**n) * Nd * dD).sum(axis=1)

def RainRate(D, Nd, Vd, dD) :
    return 6e-4 * np.pi * nthMoment(D, Vd*Nd, dD, 3)

def Dm(D,ND,dD):
  """volume weighted diameter"""
  return nthMoment(D,ND,dD,4) / nthMoment(D,ND,dD,3)

def get_exp_rr(N0, Lambda, D, v, dD, id):
  """Returns rr from a file or computes it again"""
  if os.path.isfile(id):
    rr = io.loadmat(id)['rr']
  else:
    N0 = (np.tile(N0, (100,1))).T
    Lambda = (np.tile(Lambda, (100,1))).T
    ND = N0 * np.exp(-Lambda * D)
    rr = np.reshape(RainRate(D, ND, v, dD), (359,1250))
    io.savemat(id, {'rr' : rr})
  return rr


# get az, zh, abd zdr
file     = "KOUN_20050513-083020.mat"
data     = io.loadmat(file)
data     = data['RawData'][0,0]
Zh_cond  = (data['Zh']== -99900.0) | (data['Zh'] == -32768.0)
Zdr_cond = (data['Zdr']== -99900) | (data['Zdr'] == -32768) | (data['Zdr'] > 4.7) | (data['Zdr'] < -4.7)
Zh       = np.ma.masked_array(data['Zh'], Zh_cond )
Zdr      = np.ma.masked_array(data['Zdr'], Zdr_cond)
data     = np.load('fwscamp.npz')
fa_pi    = data['fapi']
fb_pi    = data['fbpi']
fa_0     = data['fa0']
fb_0     = data['fb0']

# this is the classifcation result from the prev. project
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
Zh       = 10**(Zh/10)
Zdr      = np.ma.masked_array(Zdr,mask)

D,dD     = np.linspace(0.08,8,100,retstep=True)
caplam   = np.linspace(0,20,100)
caplam   = np.tile(caplam, (caplam.shape[0],1))
# transpose so the correct axis is used in future computations
caplam   = np.transpose(caplam)
caplam1d = caplam[:,0]

#exp model first. We omit N0 because it cancels out in get_Zdr_lambda
ND =  np.exp(-caplam * D)
Zdr_lambda = get_Zdr_lambda(ND, dD, fa_pi, fb_pi)
#reverse the function to be increasing so we can interpolate
sort_order = np.argsort(Zdr_lambda)
fit = ip.InterpolatedUnivariateSpline(Zdr_lambda[sort_order], caplam1d[sort_order])
fitted_Zdr = fit(np.abs(Zdr.flatten()))[::-1]
#then get Lambda and remove values below -10, since we shouldn't see any that are meaningful.
exp_Lambda = np.ma.masked_array(fitted_Zdr, fitted_Zdr <= -100) 

#tile and transpose lamba so we can get one ND for each range gate/Lambda that we have
#exp_Lambda = (np.tile(exp_Lambda, (100,1))).T

#ND = np.exp(-exp_Lambda * D)
#N0 = get_N0(Zh.flatten(), Kw, lam, fa_pi, ND, dD)
#mask = (N0 <= 0) | (N0 >= 2e4)
#exp_N0 = np.ma.masked_array(N0, mask)

#plot(caplam1d, Zdr_lambda,"$Z_{dr}$ - $\Lambda$ relation, Exp model",ylabel='$Z_{dr}$ (dB)',xlabel="$\Lambda$ ($mm^{-1}$)",figname="zdr-lam-exp.png")
#hist(exp_Lambda.compressed(),1000, "$\Lambda$ distribution for given data, Exp. model", "$\Lambda$", "Number", "all-caplam.png", axis=[-5,15,0,1.8e6])
#hist(exp_N0.compressed(),300, "$N_0$ distribution for given data, Exp. model", "$N_0$", "Number", "N0.png")

# now C-G ND. Except for another parameter mu, it's the same as above.
mu = -0.0201 * caplam1d**2 + 0.902*caplam1d - 1.718
ND = D**mu * np.exp(-caplam * D)
Zdr_lambda = get_Zdr_lambda(ND, dD, fa_pi, fb_pi)
sort_order = np.argsort(Zdr_lambda)
fit = ip.InterpolatedUnivariateSpline(Zdr_lambda[sort_order], caplam1d[sort_order])
fitted_Zdr = fit(np.abs(Zdr.flatten()))[::-1]
cg_Lambda = np.ma.masked_array(fitted_Zdr, fitted_Zdr <= -10) 


#cg_Lambda = np.tile(cg_Lambda, (100,1))
#mu = -0.0201 * cg_Lambda**2 + 0.902*cg_Lambda - 1.718
#a = D**mu
#ND = a * np.exp(-cg_Lambda * D)
#N0 = get_N0(Zh.flatten(), Kw, lam, fa_pi, ND, dD)
#cg_N0 = np.ma.masked_array(N0, N0 <= -100)

#plot(caplam1d, Zdr_lambda, "$Z_{dr}$ - $\Lambda$ relation, C-G model",ylabel='$Z_{dr}$ (dB)', xlabel="$\Lambda$ ($mm^{-1}$)",figname="zdr-lam-cg.png")
#hist(cg_Lambda.compressed(),1000, "$\Lambda$ distribution for given data, C-G model", "$\Lambda$", "Number", "all-caplam-cg.png", axis=[1,10,0,8000])
#hist(cg_N0.compressed(),300, "$N_0$ distribution for given data, C-G model", "$N_0$", "Number", "N0-cg.png")


#sys.exit()
# Rainfall rate and D_m
v = -0.1021 + 4.932*D- 0.9551*D**2 + 0.07934*D**3 - 0.002362*D**4

# get rainfall rate for exponential ND
# that means we need to tile N0 and Lambda so we end up with 
# gates x 100 matrix: A 100 data point ND for each range gate
#print exp_N0.max()
#print exp_N0.min()
#print exp_Lambda.max()
#print exp_Lambda.min()
#exp_rr = get_exp_rr(exp_N0, exp_Lambda, D, v, dD, 'cache-exp-rr.mat')



## convert our polar info into cartesian
#el = 0.5
#el_rad = el/180*np.pi

## each range gate is 250m
#del_r = .250
#xx = np.arange(0, exp_rr.shape[1]) * del_r
#az = np.arange(359)
#yy = az/180 * np.pi
#r,az_rad = np.meshgrid(xx,yy)

#x = r*np.cos(el_rad) * np.sin(az_rad)
#y = r*np.cos(el_rad) * np.cos(az_rad)
#z = r*np.sin(el_rad)

#fig = plt.figure(figsize=(15,9));
#ax = plt.axes()
#cax = ax.pcolor(exp_rr)
#cb = fig.colorbar(cax)
##cb.ax.set_yticklabels(["Nothing", "GC", "BS", "DS", "WS", "CR", "GR", "BD", "RA", "HR", "RH"])
#ax.axis([-150, 150, -150, 150])
#ax.set_title("Rain rate, El=$0.5^{\circ}$")
#ax.set_ylabel(r'meridonal distance/km')
#ax.set_xlabel("zonal distance/km")
#plt.savefig("exp-rr.png")

#mu   = np.tile(mu, (all_caplam.shape[0],1))

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
koun_t = np.array(np.round((koun_t - 7) * 60), dtype='int');
# convert to linear
koun_Zh = 10**(np.array(l_Zh, dtype=float)/10)
koun_Zdr = np.array(l_Zdr, dtype=float)
# just use some single example values
#koun_Zh = 10**(40/10)
#koun_Zdr = 2

D,dD     = np.linspace(0.08,8,100,retstep=True)

# figure out Lambda for each given data point from KOUN
fitted_Zdr = fit(np.abs(koun_Zdr))
#fitted_Zdr = fit(np.abs(koun_Zdr))[::-1]
Lambda1d = np.ma.masked_array(fitted_Zdr, fitted_Zdr <= -10) 
Lambda = np.ma.masked_array(fitted_Zdr, fitted_Zdr <= -10) 
mu1d = -0.0201 * Lambda**2 + 0.902*Lambda - 1.718
Lambda = (np.tile(Lambda, (100,1))).T
print Lambda
#print D.shape

# do the same for N0 and mask out very negative values
ND = np.exp(-Lambda * D)
N0 = get_N0(koun_Zh, Kw, lam, fa_pi, ND, dD)
N0 = np.ma.masked_array(N0, N0 <= -100)

print N0


# These data are a little too sparse to do a meaningful histogram of, so just give a list
#print "Lambda values for KOUN data"
#print Lambda
#print "N_0 values for KOUN data"
#print N0

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

# get requested products for the DSD measured data
dD_meas = 0.2
v_meas = -0.1021 + 4.932*D_meas - 0.9551*D_meas**2 + 0.07934*D_meas**3 - 0.002362*D_meas**4
R_calc = RainRate(D_meas, ND_meas, v_meas, dD_meas)
Dt_meas = nthMoment(D_meas, ND_meas, dD_meas, 0)
Dm_meas = Dm(D_meas, ND_meas, dD)

# repeat for the retrieved data from the radar
#a = (N0 * (D**mu).T).T
ND = N0 * np.exp(-Lambda * D)
v = -0.1021 + 4.932*D- 0.9551*D**2 + 0.07934*D**3 - 0.002362*D**4
R_retr = RainRate(D, ND, v, dD)
Dt_retr = nthMoment(D, ND, dD, 0)
Dm_retr = Dm(D, ND, dD)


# plot side by side
fig = plt.figure(figsize=(15,9));
ax = plt.axes()
ax.plot(Dt_meas[koun_t])
ax.plot(Dt_retr)
ax.legend(["meas", 'retr'])
ax.set_title("Total Number Concentration")
ax.set_xlabel("t")
ax.set_ylabel("$N_t$ (#)")
plt.savefig("2-Dt.png")

fig = plt.figure(figsize=(15,9));
ax = plt.axes()
ax.plot(R_calc[koun_t])
ax.plot(R_retr)
ax.legend(["meas", 'retr'])
ax.set_title("Rainfall rate")
ax.set_xlabel("t")
ax.set_ylabel("R ($mm*h^-1$)")
plt.savefig("2-R.png")

fig = plt.figure(figsize=(15,9));
ax = plt.axes()
ax.plot(ND_meas[koun_t])
#ax.plot(ND)
ax.legend(["meas", 'retr'])
ax.set_title("N(D)")
ax.set_xlabel("t")
ax.set_ylabel("")
plt.savefig("2-ND.png")

fig = plt.figure(figsize=(15,9));
ax = plt.axes()
ax.plot(ND)
#ax.plot(ND)
ax.legend(["meas", 'retr'])
ax.set_title("N(D)")
ax.set_xlabel("t")
ax.set_ylabel("")
plt.savefig("2-ND_ret.png")

fig = plt.figure(figsize=(15,9));
ax = plt.axes()
ax.plot(Dm_meas[koun_t])
ax.plot(Dm_retr)
ax.plot((mu1d+4)/Lambda1d)
ax.legend(["meas", 'retr', '$\mu$ * 4 / lambda method'])
ax.set_title("mean volume diameter")
ax.set_xlabel("t")
ax.set_ylabel("$D_m$ (mm)")
plt.savefig("2-Dm.png")
