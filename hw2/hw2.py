import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as consts

###############################################
############## Needed Functions ###############
###############################################

def e_inf(water_or_ice, temp):
  """Calculate epsilon_infinity for water or ice, depending on temperature"""
  if water_or_ice == "water":
    return 5.27137 + 0.021647*temp - 0.00131198*temp**2
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
    return 0.0009990288*np.exp(13200/((temp+273)*1.9869))

def sigma(water_or_ice, temp):
  """Calculate sigma"""
  if water_or_ice == "water":
    return 12.5664e8
  else:
    return 1.26 * np.exp(-12500/((temp+273)* 1.9869))

def e_real(water_or_ice, lam, temp):
  """Calculate the real part of the dielectric constant"""
  lambda_s_local = lambda_s(water_or_ice, temp)
  alpha_local = alpha(water_or_ice, temp)
  e_inf_local = e_inf(water_or_ice, temp)
  return e_inf_local + ((e_s(water_or_ice, temp) - e_inf_local) * (1 + (lambda_s_local/lam)**(1-alpha_local) * np.sin(alpha_local*np.pi/2))) / \
      (1 + 2*(lambda_s_local/lam)**(1-alpha_local) * np.sin(alpha_local*np.pi/2) + (lambda_s_local/lam)**(2-2*alpha_local))

def e_img(water_or_ice, lam, temp):
  """Calculate the real part of the dielectric constant"""
  lambda_s_local = lambda_s(water_or_ice, temp)
  alpha_local = alpha(water_or_ice, temp)
  e_inf_local = e_inf(water_or_ice, temp)
  return (((e_s(water_or_ice, temp) - e_inf_local)*(lambda_s_local/lam)**(1-alpha_local) * np.cos(alpha_local*np.pi/2))/ \
      (1 + 2*(lambda_s_local/lam)**(1-alpha_local) * np.sin(alpha_local * np.pi/2) + (lambda_s_local/lam)**(2-2*alpha_local))) + sigma(water_or_ice,temp)*lam / (18.8496e10)


# generated this formula using sympy:
# solve( , x)
def PS_mixing(e1,e2, f1,f2):
  """Calculate the effective dielectric constant for two media mixed together using Polder-van Sanden"""

  # solve the Polder-van Sanden mixing formula by hand for the coefficients
  # (f1*e1 - f1*x)/(e1 + 2*x) + (f2*e2 - f2*x)/(e2+2*x) = 0
  a = (-2*f1-2*f2)
  b = 2*e1*f1 + 2*e2*f2 - f1*e2 - e1*f2
  c = e1*e2*f1+ e1*e2*f2

  t = np.sqrt(b**2 - 4*a*c)
  return (-b + t)/(2*a), (-b - t)/(2*a)

def MG_mixing(e1,e2,f):
  """Maxwell-Garnet mixing formula"""
  y = (e2-e1)/(e2+2*e1)
  return (1+2*f*y)*e1/(1-f*y)

############################################################### 
######################## MAIN function ########################
############################################################### 

if __name__ == '__main__' :
  
  ############## One ################

  # wavelength at S band
  lS = consts.c / (2.8e9) * 100
  # wavelength at C band
  lC = consts.c / (5.6e9) * 100
  # wavelength at X band
  lX = consts.c / (10e9) * 100
  # wavelength at Ka band
  lKa = consts.c / (35e9) * 100


  t = (0, 10, 20)
  l = (lS, lC, lX, lKa)
  labels = ("S-band ", "C-band ", "X-band ", "Ka-band")

  print("Water")
  print

  for theT in t:
    for i in range(4):
      result = np.complex(e_real("water", l[i], theT), e_img("water", l[i], theT))
      print "%s @ %2.0f degrees: %s" % (labels[i],theT,result)

  print

  ################# TWO ##########################


  t = (-20, -10, 0)

  print("Ice")
  print

  for theT in t:
    for i in range(4):
      result = np.complex(e_real("ice", l[i], theT), e_img("ice", l[i], theT))
      print "%s @ %3.0f degrees: %s" % (labels[i],theT,result)



  ################# THREE ####################

  # snow density
  rho_s = 0.2

  # ice density
  rho_i = 0.9167

  # water density
  rho_w = 1

  # first dry snow, 0 melting rate
  gamma_w = 0

  # find fractional volumes 
  f_w = gamma_w * rho_s / rho_w
  f_i = rho_s / rho_i
  f_a = 1 - f_i - f_w

  # get the dielectric constant for ice and water at 0 deg, S-band again
  e_i = np.complex(e_real("ice", lS, 0), e_img("ice", lS, 0))
  e_w = np.complex(e_real("water", lS, 0), e_img("water", lS, 0))

  # dielectric constant for dry snow
  e_s = PS_mixing(1.0,e_i,f_a,f_i)

  # only the second solution makes sense
  e_s = e_s[1]
  print
  print "e_s = %s" % e_s
  print
  # melting rate from 0 to 100%
  gamma_w = np.linspace(0,0.5,100)

  # find fractional volumes.
  # here we have water and snow, so the fractional volume of snow
  # is 1 - whatever part water takes up. We don't care what 
  # the composition of snow is in this case.
  f_w = gamma_w * rho_s / rho_w
  f_s = 1 - f_w

  # compute the dielectric constant for melting snow, snow background
  # and water inclusion
  e_m1 = MG_mixing(e_w, e_s, f_s)
  
  # compute the dielectric constant for melting snow, water background
  # and snow inclusion
  e_m2 = MG_mixing(e_s, e_w, f_w)

  # for good measure, calculate the PS mixing ration as well
  e_m3 = PS_mixing(e_w, e_s, f_w, f_s)


  # plot both real and imaginary parts of the result
  fig = plt.figure(figsize=(15,9));
  ax = fig.add_subplot(2,1,1)
  ax.plot(gamma_w * 100, np.real(e_m1), gamma_w * 100, np.real(e_m2), gamma_w * 100, np.real(e_m3[1]))
  ax.legend(["MG: water/snow", "MG: snow/water", "PS: snow + water"], loc="upper left")
  ax.set_ylabel("real($\epsilon_e$)")
  ax.set_title("Dielectric constant, real(top) and imaginary(bottom)")

  ax = fig.add_subplot(2,1,2)
  ax.plot(gamma_w * 100, np.imag(e_m1), gamma_w * 100, np.imag(e_m2), gamma_w * 100, np.imag(e_m3[1]))
  ax.legend(["MG: water/snow", "MG: snow/water", "PS: snow + water"], loc="upper left")
  ax.set_xlabel("% melted")
  ax.set_ylabel("imag($\epsilon_e$)")

  plt.savefig("P3.png")
  plt.show()






