#!/usr/bin/env python

import numpy as np
import sys
from numpy import tile
import matplotlib as m
import matplotlib.pyplot as plt
import scipy.io as io
import peach.fuzzy.mf as mf


def f1(Z):
  """compute f1 to be used in the MF"""
  return -0.5 + 2.5e-3 * Z + 7.5e-4 * Z**2

def f2(Z):
  """compute f2 to be used in the MF"""
  return 0.68 - 4.81e-2 * Z + 2.92e-3 * Z**2

def f3(Z):
  """compute f3 to be used in the MF"""
  return 1.42 + 6.67e-2 * Z + 4.85e-4 * Z**2

def calcFuzzies():
  """compute fuzzy values for each of the 10 categories"""
  file       = "KOUN_20050513-083020.mat"
  data       = io.loadmat(file)

  data       = data['RawData'][0,0]

  Az         = data['Az']
  Zh         = np.ma.masked_array(data['Zh'], data['Zh'] < 0)
  Zdr        = np.ma.masked_array(data['Zdr'], data['Zdr'] < 0)
  rhoHV      = np.ma.masked_array(data['RhoHV'], data['RhoHV'] < 0)
  #phiDP     = data['PhiDP']
  #Kdp       = data['Kdp']
  HydroClass = data['HydroClass']
  time       = data['time']
  gatewidth  = data['Gatewidth']
  beamwidth  = data['Beamwidth']

  #print Az.shape
  #print Zh.shape
  #print Zdr.shape
  #print rhoHV.shape
  #print phiDP.shape
  #print Kdp.shape
  #print HydroClass.shape
  #print time.shape
  #print gatewidth.shape
  #print beamwidth.shape

  f1 = f1(Zh)
  f2 = f2(Zh)
  f3 = f3(Zh)

  GC_Z = mf.Trapezoid(15, 20, 70, 80)(Zh)
  BS_Z = mf.Trapezoid(5, 10, 20, 30)(Zh)
  DS_Z = mf.Trapezoid(5, 10, 35, 40)(Zh)
  WS_Z = mf.Trapezoid(25, 30, 40, 50)(Zh)
  CR_Z = mf.Trapezoid(0, 5, 20, 25)(Zh)
  GR_Z = mf.Trapezoid(25, 35, 50, 55)(Zh)
  BD_Z = mf.Trapezoid(20, 25, 45, 50)(Zh)
  RA_Z = mf.Trapezoid(5, 10, 45, 50)(Zh)
  HR_Z = mf.Trapezoid(40, 45, 55, 60)(Zh)
  RH_Z = mf.Trapezoid(45, 50, 75, 80)(Zh)


  GC_Zdr = mf.Trapezoid(-4, -2, 1, 2)(Zdr)
  BS_Zdr = mf.Trapezoid(0, 2, 10, 12)(Zdr)
  DS_Zdr = mf.Trapezoid(-0.3, 0, 0.3, 0.6)(Zdr)
  WS_Zdr = mf.Trapezoid(0.5, 1, 2, 3)(Zdr)
  CR_Zdr = mf.Trapezoid(0.1, 0.4, 3.0, 3.3)(Zdr)


  GC_rho = mf.Trapezoid(0.5, 0.6, 0.9, 0.95)(rhoHV)
  BS_rho = mf.Trapezoid(0.3, 0.5, 0.8, 0.83)(rhoHV)
  DS_rho = mf.Trapezoid(0.95, 0.98, 1, 1.01)(rhoHV)
  WS_rho = mf.Trapezoid(0.88, 0.92, 0.95, 0.985)(rhoHV)
  CR_rho = mf.Trapezoid(0.95, 0.98, 1, 1.01)(rhoHV)
  GR_rho = mf.Trapezoid(0.9, 0.97, 1, 1.01)(rhoHV)
  BD_rho = mf.Trapezoid(0.92, 0.95, 1, 1.01)(rhoHV)
  RA_rho = mf.Trapezoid(0.95, 0.97, 1, 1.01)(rhoHV)
  HR_rho = mf.Trapezoid(0.92, 0.95, 1, 1.01)(rhoHV)
  RH_rho = mf.Trapezoid(0.85, 0.9, 1, 1.01)(rhoHV)

  GR_Zdr = np.zeros((Zh.shape[0], Zh.shape[1]))
  BD_Zdr = np.zeros((Zh.shape[0], Zh.shape[1]))
  RA_Zdr = np.zeros((Zh.shape[0], Zh.shape[1]))
  HR_Zdr = np.zeros((Zh.shape[0], Zh.shape[1]))
  RH_Zdr = np.zeros((Zh.shape[0], Zh.shape[1]))


  for i in range(Zh.shape[0]) :
    sys.stdout.write("\r") 
    sys.stdout.write(str(i)) 
    sys.stdout.flush()
    for j in range(Zh.shape[1]):

      z = Zh[i][j]
      GR_Zdr[i][j] = mf.Trapezoid(-0.3, 0, f1[i][j], f1[i][j] + 0.3)(Zdr[i][j])
      BD_Zdr[i][j] = mf.Trapezoid(f2[i][j] - 0.3, f2[i][j], f3[i][j], f3[i][j]+1.)(Zdr[i][j])
      RA_Zdr[i][j] = mf.Trapezoid(f1[i][j]-0.3, f1[i][j], f2[i][j], f2[i][j]+0.5)(Zdr[i][j])
      HR_Zdr[i][j] = mf.Trapezoid(f1[i][j]-0.3, f1[i][j], f2[i][j], f2[i][j]+0.5)(Zdr[i][j])
      RH_Zdr[i][j] = mf.Trapezoid(-0.3, 0, f1[i][j], f1[i][j]+0.5)(Zdr[i][j])

  d = {}
  d['GC'] = GC_Z + GC_Zdr + GC_rho
  d['BS'] = BS_Z + BS_Zdr + BS_rho
  d['DS'] = DS_Z + DS_Zdr + DS_rho
  d['WS'] = WS_Z + WS_Zdr + WS_rho
  d['CR'] = CR_Z + CR_Zdr + CR_rho
  d['GR'] = GR_Z + GR_Zdr + GR_rho
  d['BD'] = BD_Z + BD_Zdr + BD_rho
  d['RA'] = RA_Z + RA_Zdr + RA_rho
  d['HR'] = HR_Z + HR_Zdr + HR_rho
  d['RH'] = RH_Z + RH_Zdr + RH_rho

  io.savemat('res', d)



# get az, zh, abd zdr
file       = "KOUN_20050513-083020.mat"
data       = io.loadmat(file)
data       = data['RawData'][0,0]
az         = np.transpose(data['Az'])[0]
Zh         = np.ma.masked_array(data['Zh'], data['Zh'] < 0)
Zdr        = np.ma.masked_array(data['Zdr'], data['Zdr'] < 0)

data = io.loadmat('res.mat')
dim0 = data['GC'].shape[0]
dim1 = data['GC'].shape[1]

# 5 regions to compute for
idx0 = 110000 / 250
idx1 = 151000 / 250
idx2 = 180000 / 250
idx3 = 238000 / 250
idx4 = dim1



# first, 0-110km, 6 products
first = np.zeros((11, dim0, dim1))
first[0] = np.zeros((dim0, dim1))
first[1] = data['GC']
first[2] = data['BS']
first[3] = np.zeros((dim0,dim1))
first[4] = np.zeros((dim0,dim1))
first[5] = np.zeros((dim0,dim1))
first[6] = np.zeros((dim0,dim1))
first[7] = data['BD']
first[8] = data['RA']
first[9] = data['HR']
first[10] = data['RH']

first = first[:,:,0:idx0]
first_res = np.argmax(first, axis=0)

# second, 110km-151km, 8 products
second = np.zeros((11, dim0, dim1))
second[0] = np.zeros((dim0, dim1))
second[1] = data['GC']
second[2] = data['BS']
second[3] = np.zeros((dim0,dim1)) # DS
second[4] = data['WS']
second[5] = np.zeros((dim0,dim1)) # CR
second[6] = data['GR']
second[7] = data['BD']
second[8] = data['RA']
second[9] = data['HR']
second[10] = data['RH']

second = second[:,:,idx0:idx1]
second_res = np.argmax(second, axis=0)

# third, 110km-151km, 7 products
third = np.zeros((11, dim0, dim1))
third[0] = np.zeros((dim0, dim1))
third[1] = data['GC']
third[2] = data['BS']
third[3] = data['DS'] 
third[4] = data['WS']
third[5] = np.zeros((dim0,dim1)) # CR
third[6] = data['GR']
third[7] = data['BD']
third[8] = np.zeros((dim0, dim1)) # RA
third[9] = np.zeros((dim0, dim1)) # HR
third[10] = data['RH']

third = third[:,:,idx1:idx2]
third_res = np.argmax(third, axis=0)

# fourth, 110km-151km, 8 products
fourth = np.zeros((11, dim0, dim1))
fourth[0] = np.zeros((dim0, dim1))
fourth[1] = data['GC']
fourth[2] = data['BS']
fourth[3] = data['DS'] 
fourth[4] = data['WS']
fourth[5] = data['CR']
fourth[6] = data['GR']
fourth[7] = data['BD']
fourth[8] = np.zeros((dim0, dim1)) # RA
fourth[9] = np.zeros((dim0, dim1)) # HR
fourth[10] = data['RH']

fourth = fourth[:,:,idx2:idx3]
fourth_res = np.argmax(fourth, axis=0)


# fifth, 110km-151km, 4 products
fifth = np.zeros((11, dim0, dim1))
fifth[0] = np.zeros((dim0, dim1))
fifth[1] = np.zeros((dim0, dim1)) # GC
fifth[2] = np.zeros((dim0, dim1)) # BS
fifth[3] = data['DS'] 
fifth[4] = np.zeros((dim0, dim1)) # WS
fifth[5] = data['CR']
fifth[6] = data['GR']
fifth[7] = np.zeros((dim0, dim1)) # BD
fifth[8] = np.zeros((dim0, dim1)) # RA
fifth[9] = np.zeros((dim0, dim1)) # HR
fifth[10] = data['RH']

fifth = fifth[:,:,idx3:idx4]
fifth_res = np.argmax(fifth, axis=0)


total = np.concatenate((first_res, second_res, third_res, fourth_res, fifth_res), axis=1)

# convert our polar info into cartesian
el = 0.5
el_rad = el/180*np.pi

# each range gate is 250m
del_r = .250
xx = np.arange(0, dim1) * del_r
yy = az/180 * np.pi
r,az_rad = np.meshgrid(xx,yy)

x = r*np.cos(el_rad) * np.sin(az_rad)
y = r*np.cos(el_rad) * np.cos(az_rad)
z = r*np.sin(el_rad)

# define colormap
colors = np.zeros((11,3))
colors[0] = [1,1,1]
colors[1] = [93./255,1,1]
colors[2] = [15./255,0,204./255]
colors[3] = [58./255,207./255,0]
colors[4] = [42./255,155./255,0]
colors[5] = [1,1,50./255]
colors[6] = [246./255,0,25./255]
colors[7] = [197./255,0,18./255]
colors[8] = [246./255,0,206./255]
colors[9] = [149./255,34./255,205./255]
colors[10] = [123./255,0,128./255]

fig = plt.figure(figsize=(15,9));
ax = plt.axes()
cm = m.colors.ListedColormap(list(colors))
cax = ax.pcolor(x,y,total, cmap=cm)
cb = fig.colorbar(cax)
cb.ax.set_yticklabels(["Nothing", "GC", "BS", "DS", "WS", "CR", "GR", "BD", "RA", "HR", "RH"])
ax.axis([-150, 150, -150, 150])
ax.set_title("Hydrometeor classification, El=$0.5^{\circ}$")
ax.set_ylabel(r'meridonal distance/km')
ax.set_xlabel("zonal distance/km")
plt.savefig("classification.png")

fig = plt.figure(figsize=(15,9));
ax = plt.axes()
cax = ax.pcolor(x,y,Zh)
cb = fig.colorbar(cax)
ax.axis([-150, 150, -150, 150])
ax.set_title("$Z_H \ (dBZ)$, El=$0.5^{\circ}$")
ax.set_ylabel(r'meridonal distance/km')
ax.set_xlabel("zonal distance/km")
plt.savefig("zh.png")

fig = plt.figure(figsize=(15,9));
ax = plt.axes()
cax = ax.pcolor(x,y,Zdr)
cb = fig.colorbar(cax)
ax.axis([-150, 150, -150, 150])
ax.set_title("$Z_H \ (dBZ)$, El=$0.5^{\circ}$")
ax.set_ylabel(r'meridonal distance/km')
ax.set_xlabel("zonal distance/km")
plt.savefig("zdr.png")
plt.show()
