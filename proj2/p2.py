#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.io as io


file = "KOUN_20050513-083020.mat"
data = io.loadmat(file)

data = data['RawData'][0,0]

Az = data['Az']
Zh = data['Zh']
Zdr = data['Zdr']
rhoHV = data['RhoHV']
phiDP = data['PhiDP']
Kdp = data['Kdp']
HydroClass = data['HydroClass']
time = data['time']
gatewidth = data['Gatewidth']
beamwidth = data['Beamwidth']

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
