#!/usr/bin/env python

import numpy as np
import sys
from numpy import tile
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




data = io.loadmat('res.mat')
#fig = plt.figure(figsize=(15,9));
#ax = plt.axes()
#ax.plot(foo)
#ax.plot(testdata)
##ax.set_title("$K_{dp}$")
##ax.set_ylabel(r'$deg/km$')
##ax.set_xlabel("rainfall rate R (mm/hr)")
#ax.legend(['fuzzy', 'Zh'])
##ax.set_yscale('log')
#plt.savefig("test.png")
##plt.show()
