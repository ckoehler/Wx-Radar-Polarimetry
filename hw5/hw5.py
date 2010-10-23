#!/usr/bin/env python

import numpy as np
import scipy.constants as consts
import matplotlib.pyplot as plt
import scipy.io as io
import scipy.stats as stats


# load data from the simulation and assign to local variables
data = io.loadmat('data.mat')
A = data['A'][0]
Ib = data['Ib'][0]
I = data['I'][0]
Q = data['Q'][0]
phi = data['phi'][0]

#######################
##### P3 check A ######
#######################
plt.figure(figsize=(15,9));
ax = plt.axes()

# get rayleigh fit
fit_args = stats.rayleigh.fit(A)

# plot histogram and normalize it
n, bins, patches = ax.hist(A, 100, normed=True)

# compute the range of the rayleigh pdf
theRange = np.linspace(bins.min(), bins.max(), np.ceil(bins.max()))

# figure out the pdf
pdf = stats.rayleigh.pdf(theRange, loc=fit_args[0], scale=fit_args[1])

# plot
ax.plot(theRange, pdf)
ax.set_title("Histogram of $A$ with fitted Rayleigh PDF")
ax.set_xlabel("Amplitudes $A$")
ax.set_ylabel("# of Occurances, normalized")
plt.savefig("p3-A.png")

#########################
##### P3 check I_b ######
#########################
plt.figure(figsize=(15,9));
ax = plt.axes()

# get rayleigh fit
fit_args = stats.expon.fit(Ib)

# plot histogram and normalize it
n, bins, patches = ax.hist(Ib, 100, normed=True)

# compute the range of the rayleigh pdf
theRange = np.linspace(bins.min(), bins.max(), np.ceil(bins.max()))

# figure out the pdf
pdf = stats.expon.pdf(theRange, loc=fit_args[0], scale=fit_args[1])

# plot
ax.plot(theRange, pdf)
ax.set_title("Histogram of $I_b$ with fitted Exponential PDF")
ax.set_xlabel("Intensities $I_b$")
ax.set_ylabel("# of Occurances, normalized")
plt.savefig("p3-Ib.png")

#######################
##### P3 check I ######
#######################
plt.figure(figsize=(15,9));
ax = plt.axes()

# get rayleigh fit
fit_args = stats.norm.fit(I)

# plot histogram and normalize it
n, bins, patches = ax.hist(I, 100, normed=True)

# compute the range of the rayleigh pdf
theRange = np.linspace(bins.min(), bins.max(), np.ceil(bins.max()))

# figure out the pdf
pdf = stats.norm.pdf(theRange, loc=fit_args[0], scale=fit_args[1])

# plot
ax.plot(theRange, pdf)
ax.set_title("Histogram of $I$ with fitted Gaussian PDF")
ax.set_xlabel("In-phase values $I$")
ax.set_ylabel("# of Occurances, normalized")
plt.savefig("p3-I.png")


#######################
##### P3 check Q ######
#######################
plt.figure(figsize=(15,9));
ax = plt.axes()

# get rayleigh fit
fit_args = stats.norm.fit(Q)

# plot histogram and normalize it
n, bins, patches = ax.hist(Q, 100, normed=True)

# compute the range of the rayleigh pdf
theRange = np.linspace(bins.min(), bins.max(), np.ceil(bins.max())) 

# figure out the pdf
pdf = stats.norm.pdf(theRange, loc=fit_args[0], scale=fit_args[1])

# plot
ax.plot(theRange, pdf)
ax.set_title("Histogram of $Q$ with fitted Gaussian PDF")
ax.set_xlabel("Quadrature values $Q$")
ax.set_ylabel("# of Occurances, normalized")
plt.savefig("p3-Q.png")

#########################
##### P3 check phi ######
#########################
plt.figure(figsize=(15,9));
ax = plt.axes()

# get rayleigh fit
fit_args = stats.norm.fit(Q)

# plot histogram and normalize it
n, bins, patches = ax.hist(phi, 100, normed=True)

# compute the range of the rayleigh pdf
theRange = np.linspace(bins.min(), bins.max(), np.ceil(bins.max())) 

# figure out the pdf
pdf = stats.norm.pdf(theRange, loc=fit_args[0], scale=fit_args[1])
pdf = 1./2./np.pi * np.ones(np.ceil(bins.max())) 

# plot
ax.plot(theRange, pdf)
ax.set_title("Histogram of $\phi$ with fitted constant PDF")
ax.set_xlabel("Phase values $\phi$")
ax.set_ylabel("# of Occurances, normalized")
plt.savefig("p3-phi.png")


######################################
##### P4 check expected values  ######
######################################
sigma2 = np.mean(I**2)
sigma = np.sqrt(sigma2)
I_mean = np.mean(I)
print "<I> = %f " % I_mean

Q_mean = np.mean(Q)
print "<Q> = %f " % Q_mean

A_mean = np.mean(A)
print "<A> = %f " % A_mean
print "<A> calculated = %f " % (sigma * np.sqrt(np.pi/2))

Ib_mean = np.mean(Ib)
print "<Ib> = %f " % Ib_mean
print "<Ib> calculated = %f " % (sigma**2 * 2)

phi_mean = np.mean(phi)
print "<phi> = %f " % phi_mean

#plt.show()
