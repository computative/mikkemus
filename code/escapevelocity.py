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

h = 1e-3
N = 10000
G = 4*pi**2
M = 1
x = []
y = []
vs = append(linspace(0,7.5,5),linspace(7,11,20))
vs = append(vs,linspace(11,20,5))
for v0 in vs:#linspace(pi,4*pi,10):
    r = zeros((N,2))
    v = zeros((N,2))
    r[0] = [1,0]
    v[0] = v0*array([1, 0*0.01])
    print v0
    
    # velocity verlet integration
    for i in range(N-1):
        r[i+1] = r[i] + h*v[i] - (h**2/2)*G*M*r[i]/norm(r[i])**3
        v[i+1] = v[i] - G*M*(h/2)*(r[i+1]/norm(r[i+1])**3 + r[i]/norm(r[i])**3 )
        if (v[i+1,0] < 0): break
    y.append((v[0,0]**2 - v[argmin(v[:,0])-1,0]**2)**.5 )
    x.append(v0)


import matplotlib.pyplot as plt

f, ax = plt.subplots(figsize=(0.8*4,0.8*3))

ax.plot(x,y, 'k');
ax.plot([0,20],[8**0.5*pi,8**.5*pi], 'k:');
ax.set_title(ur'Escape velocity for planet',fontsize=11);
ax.set_xlabel(ur'$v_0$ \quad [AU/yr]');
ax.set_ylabel(ur'$\big(v_0^2 - \|\vec{v}\|_\infty^2\big)^{1/2}$ \quad [AU/yr]');
ax.set_ylim(0,9*1.6)
ax.set_yticks([8**.5*pi])
ax.text(1,8**.5*pi+0.5,ur'$(2GM/r)^{1/2}$',fontsize=11)
ax.set_yticklabels([ur''],rotation='vertical')
ax.set_xticks([0,4,8,12,16,20])
plt.tight_layout(0.5)
plt.savefig("/home/marius/Dokumenter/fys4150/project3/escapevelocity.pgf")
plt.savefig("/home/marius/Dokumenter/fys4150/project3/escapevelocity.png")
