# -*- coding: utf-8 -*-
from numpy import *
from numpy.linalg import *
import matplotlib.pyplot as plt

# load result of GR simulation
a = loadtxt('../resources/steps.txt')
print a[:201:50]*206264.806


# plot earth sun solution
f, ax = plt.subplots(figsize=(0.8*4,0.8*3))
ax.plot([0,50,100],a[:201:100]*206264.806,  'k', linewidth=1);
ax.set_title(ur'Regression perihelion precession' ,fontsize=11);
ax.set_xlabel(ur'Time    [yrs]');
ax.set_ylabel(ur'Precession    [arcsec]'); 
plt.tight_layout(0.5)
plt.savefig("../benchmark/grcorrection.png")


# plot earth sun solution
f, ax = plt.subplots(figsize=(0.8*4,0.8*3))
ax.plot(linspace(0,100, len(a)),a*206264.806,  'k', linewidth=1);
ax.set_title(ur'Raw data, perihelion precession' ,fontsize=11);
ax.set_xlabel(ur'Time    [yrs]');
ax.set_ylabel(ur'Precession    [arcsec]'); 
plt.tight_layout(0.5)
plt.savefig("../benchmark/grcorrection2.png")
