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

f, ax = plt.subplots(figsize=(0.8*4,0.8*3))

# constants and initial conditions
h = 1e-3; N = int(1/h)
r = zeros((N,2)); v = zeros((N,2)); G = 4*pi**2;
M = 1; r[0] = [1,0]; v[0] = (G*M)**0.5*array([-r[0,1], r[0,0]])

# forward euler integration

for i in arange(N-1):
    r[i+1] = r[i] + h*v[i]
    v[i+1] = v[i] - h*G*M*r[i]/norm(r[i])**3

# velocity verlet integration
r2 = zeros((N,2)); v2 = zeros((N,2));r2[0] = [1,0];
v2[0] = (G*M)**0.5*array([-r[0,1], r[0,0]])

for i in range(N-1):
    r2[i+1] = r2[i] + h*v2[i] - (h**2/2)*G*M*r2[i]/norm(r2[i])**3
    v2[i+1] = v2[i] - G*M*(h/2)*(r2[i+1]/norm(r2[i+1])**3 + r2[i]/norm(r2[i]) )


# plot earth sun solution
ax.plot(r[0:499,0], r[0:499,1], 'k--', linewidth=2);
ax.plot(r2[500:999,0], r2[500:999,1], 'k');
ax.set_xlim(-1.6*4./3,1.6*4./3);
ax.set_ylim(-1.6,1.6);
ax.set_title(ur'Orbit of earth, $\Delta t =$'+str(h),fontsize=11);
ax.set_xlabel(ur'$x$ \quad [AU]');
ax.set_ylabel(ur'$y$ \quad [AU]'); 
plt.tight_layout(0.5)
plt.show()
plt.savefig("/home/marius/Dokumenter/fys4150/project3/earthsun.pgf")
plt.savefig("/home/marius/Dokumenter/fys4150/project3/earthsun.png")
