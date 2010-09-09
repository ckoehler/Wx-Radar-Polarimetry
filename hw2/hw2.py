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
    return 78.54*(1.0-4.579e-3*(temp-25.0)+1.19e-5*(temp-25)**2 - 2.8e-8 * (temp - 25)**3)
  else:
    return 203.168 + 2.5*temp + 0.15*temp**2

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
  lS = consts.c / (2.8e9)
  # wavelength at C band
  lC = consts.c / (5.6e9)
  # wavelength at X band
  lX = consts.c / (10e9)
  # wavelength at Ka band
  lKa = consts.c / (35e9)


  t = range(0,25,1)
  e_S = e_real("water", lS, np.arange(0,25.0,1))
  e_C = e_real("water", lC, np.arange(0,25.0,1))
  e_X = e_real("water", lX, np.arange(0,25.0,1))
  e_Ka = e_real("water", lKa, np.arange(0,25.0,1))

  fig = plt.figure()
  ax = fig.add_subplot(1,2,1)
  ax.plot(t, e_S, t, e_C, t, e_X, t, e_Ka)
  ax.legend(["S-band", "C-band","X-band","Ka-band"], loc="lower center")
  ax.set_xlabel("Temperature ($^\circ$C)")
  ax.set_ylabel("Dielectric constant (unitless)")
  ax.set_title("Dielectric constant, water, real part")


  e_S = e_img("water", lS, np.arange(0,25.0,1))
  e_C = e_img("water", lC, np.arange(0,25.0,1))
  e_X = e_img("water", lX, np.arange(0,25.0,1))
  e_Ka = e_img("water", lKa, np.arange(0,25.0,1))

  ax = fig.add_subplot(1,2,2)
  ax.plot(t, e_S, t, e_C, t, e_X, t, e_Ka)
  ax.legend(["S-band", "C-band","X-band","Ka-band"],loc="upper left")
  ax.set_xlabel("Temperature ($^\circ$C)")
  ax.set_ylabel("Dielectric constant (unitless)")
  ax.set_title("Dielectric constant, water, imaginary part")
  
  fig.savefig("P1.png")


  ################# TWO ##########################

  t = np.arange(-20.0,0,1)
  e_S  = e_real("ice", lS, np.arange(-20.0,0,1))
  e_C  = e_real("ice", lC, np.arange(-20.0,0,1))
  e_X  = e_real("ice", lX, np.arange(-20.0,0,1))
  e_Ka = e_real("ice", lKa, np.arange(-20.0,0,1))

  fig = plt.figure()
  ax = fig.add_subplot(1,2,1)
  ax.plot(t, e_S, t, e_C, t, e_X, t, e_Ka)
  ax.legend(["S-band", "C-band","X-band","Ka-band"], loc="upper right")
  ax.set_xlabel("Temperature ($^\circ$C)")
  ax.set_ylabel("Dielectric constant (unitless)")
  ax.set_title("Dielectric constant, ice, real part")


  e_S  = e_img("ice", lS, np.arange(-20.0,0,1))
  e_C  = e_img("ice", lC, np.arange(-20.0,0,1))
  e_X  = e_img("ice", lX, np.arange(-20.0,0,1))
  e_Ka = e_img("ice", lKa, np.arange(-20.0,0,1))

  ax = fig.add_subplot(1,2,2)
  ax.plot(t, e_S, t, e_C, t, e_X, t, e_Ka)
  ax.legend(["S-band", "C-band","X-band","Ka-band"],loc="center right")
  ax.set_xlabel("Temperature ($^\circ$C)")
  ax.set_ylabel("Dielectric constant (unitless)")
  ax.set_title("Dielectric constant, ice, imaginary part")

  fig.savefig("P2.png")


  ################# THREE ####################

  # snow density
  rho_s = 0.2

  # melting rate from 0 to 50%
  gamma_w = np.arange(0,50.0,1)



  plt.show()





