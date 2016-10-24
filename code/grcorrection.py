# -*- coding: utf-8 -*-
from numpy import *
from numpy.linalg import *

import matplotlib as mpl
mpl.use("pgf")
pgf_with_pdflatex = {
    "font.family": "serif",
    "font.serif": [],
    "font.size" : 11.0,
    "pgf.preamble": [
         r"\usepackage[utf8]{inputenc}",
         r"\usepackage[T1]{fontenc}",
         r"\usepackage{cmbright}",
         r"\usepackage{newtxtext}",
         r"\usepackage{bm}",
         r"\usepackage{amsmath,amsthm}"
         ]
}
mpl.rcParams.update(pgf_with_pdflatex)
# blue, violet, green, brown, red ,orange
mpl.rcParams['axes.color_cycle'] = ['b', '#800080', '#009900', '#964B00', '#d30704' ,'#ff6100']

import matplotlib.pyplot as plt


a = loadtxt('steps.txt')
print a[:201:50]*206264.806
f, ax = plt.subplots(figsize=(0.8*4,0.8*3))

# plot earth sun solution
ax.plot([0,50,100],a[:201:100]*206264.806,  'k', linewidth=1);
ax.set_title(ur'Regression perihelion precession' ,fontsize=11);
ax.set_xlabel(ur'Time \quad [yrs]');
ax.set_ylabel(ur'Precession \quad [arcsec]'); 
plt.tight_layout(0.5)
plt.show()
plt.savefig("/home/marius/Dokumenter/fys4150/project3/grcorrection.pgf")
plt.savefig("/home/marius/Dokumenter/fys4150/project3/grcorrection.png")


f, ax = plt.subplots(figsize=(0.8*4,0.8*3))

# plot earth sun solution
ax.plot(linspace(0,100, len(a)),a*206264.806,  'k', linewidth=1);
ax.set_title(ur'Raw data, perihelion precession' ,fontsize=11);
ax.set_xlabel(ur'Time \quad [yrs]');
ax.set_ylabel(ur'Precession \quad [arcsec]'); 
plt.tight_layout(0.5)
plt.show()
plt.savefig("/home/marius/Dokumenter/fys4150/project3/grcorrection2.pgf")
plt.savefig("/home/marius/Dokumenter/fys4150/project3/grcorrection2.png")
