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

f, ax = plt.subplots(1,3,figsize=(2*0.8*4,0.8*3))

distances = loadtxt("distances.txt")
static = loadtxt("3e.txt")
dynamic = loadtxt("3f.txt")
I = linspace(0,2*pi,len(dynamic))
exact = asarray([cos(I),sin(I)])

static = norm(static[:,:2], axis=1)
dynamic = norm(dynamic[:,:2], axis=1)
exact = norm(exact.T, axis=1)

n = static.size
freqs = fftfreq(n, d=1e-3)

from matplotlib.pyplot import *


ax[0].plot(freqs[:int(n/2.)],abs(fft(static-mean(static)))[:int(n/2.)]/n,'k')
ax[1].plot(freqs[:int(n/2.)],abs(fft(dynamic-mean(dynamic)))[:int(n/2.)]/n,'k')
ax[2].plot(freqs[:int(n/2.)],abs(fft(distances-mean(distances)))[:int(n/2.)]/n,'k')
ax[0].set_xlim(0,2)
ax[1].set_xlim(0,2)
ax[2].set_xlim(0,2)
ax[0].ticklabel_format(style = 'sci', useOffset=False, scilimits=(-2,2))
ax[1].ticklabel_format(style = 'sci', useOffset=False, scilimits=(-2,2))
ax[0].set_ylim(0,1.2*max(abs(fft(static-mean(static)))[:int(n/2.)]/n))
ax[1].set_ylim(0,1.2*max(abs(fft(dynamic-mean(dynamic)))[:int(n/2.)]/n))
ax[2].set_ylim(0,1.2*max(abs(fft(distances-mean(distances)))[:int(n/2.)]/n))
ax[0].set_ylabel(ur'Amplitude   [ AU ]')
ax[1].set_xlabel(ur'Frequencies   [ yr$^{-1}$ ]')
ax[0].set_title(ur"Static",fontsize=11)
ax[1].set_title(ur"Non-static",fontsize=11)
ax[2].set_title(ur"Earth-Jupiter distance",fontsize=11)

plt.tight_layout(0.5)
plt.savefig("/home/marius/Dokumenter/fys4150/project3/pathalteration.pgf")
plt.savefig("/home/marius/Dokumenter/fys4150/project3/pathalteration.png")
