#!/usr/bin/env python

import numpy as np
import scipy.constants as consts
import matplotlib.pyplot as plt
import scipy.io as io
import scipy.stats as stats
import scipy.special as spec


def avgA(order, sigma):
  """compute average amplitude of order "order". First order will be amplitude,
  Second order will be intensity"""
  return 2**(order/2) * sigma**order * spec.gamma(order/2 + 1)


data = io.loadmat('data.mat')
A = data['A'][0]
Ib = data['Ib'][0]
I = data['I'][0]
Q = data['Q'][0]
phi = data['phi'][0]
t = np.arange(100)


##### P3 check A ######
plt.figure(figsize=(15,9));
ax = plt.axes()

# get rayleigh fit
fit_args = stats.rayleigh.fit(A)

# plot histogram and normalize it
n, bins, patches = ax.hist(A, 100, normed=True)

# compute the range of the rayleigh pdf
theRange = np.linspace(bins.min(), bins.max(), 30)

# figure out the pdf
pdf = stats.rayleigh.pdf(theRange, loc=fit_args[0], scale=fit_args[1])

# plot
ax.plot(theRange, pdf)
ax.set_title("Amplitude A")

##### P3 check I_b ######
plt.figure(figsize=(15,9));
ax = plt.axes()

# get rayleigh fit
fit_args = stats.expon.fit(Ib)

# plot histogram and normalize it
n, bins, patches = ax.hist(Ib, 100, normed=True)

# compute the range of the rayleigh pdf
theRange = np.linspace(bins.min(), bins.max(), 800)

# figure out the pdf
pdf = stats.expon.pdf(theRange, loc=fit_args[0], scale=fit_args[1])

# plot
ax.plot(theRange, pdf)

##### P3 check I ######
plt.figure(figsize=(15,9));
ax = plt.axes()

# get rayleigh fit
fit_args = stats.norm.fit(I)

# plot histogram and normalize it
n, bins, patches = ax.hist(I, 100, normed=True)

# compute the range of the rayleigh pdf
theRange = np.linspace(bins.min(), bins.max(), 40)

# figure out the pdf
pdf = stats.norm.pdf(theRange, loc=fit_args[0], scale=fit_args[1])

# plot
ax.plot(theRange, pdf)


##### P3 check Q ######
plt.figure(figsize=(15,9));
ax = plt.axes()

# get rayleigh fit
fit_args = stats.norm.fit(Q)

# plot histogram and normalize it
n, bins, patches = ax.hist(Q, 100, normed=True)

# compute the range of the rayleigh pdf
theRange = np.linspace(bins.min(), bins.max(), 40) 

# figure out the pdf
pdf = stats.norm.pdf(theRange, loc=fit_args[0], scale=fit_args[1])

# plot
ax.plot(theRange, pdf)


##### P3 check phi ######
plt.figure(figsize=(15,9));
ax = plt.axes()

# get rayleigh fit
fit_args = stats.norm.fit(Q)

# plot histogram and normalize it
n, bins, patches = ax.hist(phi, 100, normed=True)

# compute the range of the rayleigh pdf
theRange = np.linspace(bins.min(), bins.max(), 40) 

# figure out the pdf
pdf = stats.norm.pdf(theRange, loc=fit_args[0], scale=fit_args[1])
pdf = 1./2./np.pi * np.ones(40) 

# plot
ax.plot(theRange, pdf)


##### P4 check <I>  ######
I_mean = I.mean()
print I_mean

Q_mean = Q.mean()
print Q_mean

A_mean = A.mean()
A_std = A.std()
print A_mean
print A_std * np.sqrt(np.pi/2.)


#plt.show()
