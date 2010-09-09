import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as consts

###############################################
############## Needed Functions ###############
###############################################

def e_inf(water_or_ice, temp):
  """Calculate epsilon_infinity for water or ice, depending on temperature"""
  if water_or_ice == "water":
    return 5.27137 + 0.0216474*temp - 0.00131198*temp**2
  else:
    return 3.168


def e_s(water_or_ice, temp):
  """Calculate epsilon_s for water or ice"""
  if water_or_ice == "water":
    return 78.54*(1.0-(4.579*np.power(10,-3))*(temp-25.0)+(1.19*np.power(10,-5))*(temp-25)**2 - (2.8*np.power(10,-8)) * (temp - 25)**3)
  else:
    return 203.168 + 2.5*t + 0.15*t**2

def alpha(water_or_ice, temp):
  """Calculate alpha for water or ice"""
  if water_or_ice == "water":
    return -16.8129/(temp+273)+0.0609265
  else:
    return 0.288 + 0.0052*temp + 0.00023*temp**2

def lambda_s(water_or_ice, temp):
  """Calculate lambda_s"""
  if water_or_ice == "water":
    return 0.00033836*np.exp(2513.98/(temp+273))
  else:
    return 0.288 + 0.0052*temp + 0.00023*temp**2


def e_real(water_or_ice, lam, temp):
  """Calculate the real part of the dielectric constant"""
  e_inf_local = e_inf(water_or_ice, temp)
  return e_inf_local + (e_s(water_or_ice, temp) - e_inf_local) / (1 + (lambda_s(water_or_ice, temp)/lam)**2)

def e_img(water_or_ice, lam, temp):
  """Calculate the real part of the dielectric constant"""
  lam_s = lambda_s(water_or_ice, temp)
  return ((e_s(water_or_ice, temp) - e_inf(water_or_ice, temp))*(lam_s/lam))/(1 + (lam_s/lam)**2)


 
######################## MAIN function ########################

if __name__ == '__main__' :
  
  ############## One ################

  # wavelength at S band
  lS = consts.c / (2.8*np.power(10,9))
  print "wavelength at S-band: %f" % lS
  print
  print "dielectric constant for water:"

  for temp in [0,10,20]:
    print "e real @ %2.0f degrees C: %f" % (temp, e_real("water", lS, temp))    
    print "e imag @ %2.0f degrees C: %f" % (temp, e_img("water", lS, temp))

  # wavelength at C band
  lC = consts.c / (5.6*np.power(10,9))
  print "wavelength at C-band: %f" % lC
  print
  print "dielectric constant for water:"

  for temp in [0,10,20]:
    print "e real @ %2.0f degrees C: %f" % (temp, e_real("water", lC, temp))    
    print "e imag @ %2.0f degrees C: %f" % (temp, e_img("water", lC, temp))

  # wavelength at X band
  lX = consts.c / (10*np.power(10,9))
  print "wavelength at X-band: %f" % lX
  print
  print "dielectric constant for water:"

  for temp in [0,10,20]:
    print "e real @ %2.0f degrees C: %f" % (temp, e_real("water", lX, temp))    
    print "e imag @ %2.0f degrees C: %f" % (temp, e_img("water", lX, temp))

  # wavelength at Ka band
  lKa = consts.c / (35*np.power(10,9))
  print "wavelength at Ka-band: %f" % lKa
  print
  print "dielectric constant for water:"

  for temp in [0,10,20]:
    print "e real @ %2.0f degrees C: %f" % (temp, e_real("water", lKa, temp))
    print "e imag @ %2.0f degrees C: %f" % (temp, e_img("water", lKa, temp))



