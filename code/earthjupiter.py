# -*- coding: utf-8 -*-
from numpy import *
from numpy.linalg import *
from numpy.fft import *

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

from matplotlib.pyplot import *

r = loadtxt("3e.txt")
s = loadtxt("3e-jupiter.txt")

N = len(r)


ax.plot(r[:,0],r[:,1],'k')
ax.plot(s[:,0],s[:,1],'k')
ax.scatter([0],[0],c='k')
ax.scatter(r[N-1,0],r[N-1,1],c='k')
ax.scatter(s[N-1,0],s[N-1,1],c='k')
ax.set_xlim(-4/3.*6,4/3.*6)
ax.set_ylim(-6,6)
ax.set_xlabel(ur'$x$   [ AU ]')
ax.set_ylabel(ur'$y$   [ AU ]')
ax.set_title(ur"Three-body problem, $\Delta t=0.0001$",fontsize=11)

plt.tight_layout(0.5)
plt.savefig("/home/marius/Dokumenter/fys4150/project3/earthjupiter.pgf")
plt.savefig("/home/marius/Dokumenter/fys4150/project3/earthjupiter.png")
