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

import matplotlib.pyplot as plt


# forward euler integration
for yr,j in zip([1,100],[-5,-3]):
    print "j,yr: ",j,yr
    f, ax = plt.subplots(figsize=(0.8*4,0.8*3))
    H = logspace(j,0,40)
    data = []; data2 = []
    for h in H:
        N = int(ceil(yr/h))
        print N,h
        r = zeros((N,2)); v = zeros((N,2)); G = 4*pi**2;
        M = 1; r[0] = [1,0]; v[0] = (G*M)**0.5*array([-r[0,1], r[0,0]])
        for i in arange(N-1):
            r[i+1] = r[i] + h*v[i]
            v[i+1] = v[i] - h*G*M*r[i]/norm(r[i])**3
    
        # velocity verlet integration
        r2 = zeros((N,2)); v2 = zeros((N,2));r2[0] = [1,0];
        v2[0] = (G*M)**0.5*array([-r[0,1], r[0,0]])

        for i in range(N-1):
            r2[i+1] = r2[i] + h*v2[i] - (h**2/2)*G*M*r2[i]/norm(r2[i])**3
            v2[i+1] = v2[i] - G*M*(h/2)*(r2[i+1]/norm(r2[i+1])**3 + r2[i]/norm(r2[i])**3 )
        data.append( max( norm(r, axis=1)-1 ) )
        data2.append( max( norm(r2, axis=1)-1 ) )

    ax.loglog(H,data,'k', label='Forward Euler')
    ax.loglog(H,data2,'k--', label='Velocity Verlet')
    ax.set_xlabel(ur'$\Delta t$ \quad [yrs]')
    ax.set_ylabel(ur'Error $\|\vec{r} - \vec{s}\|_\infty$ \quad [AU]')
    ax.set_title('Error bound after '+str(yr)+' yrs simulation', fontsize=11)
    plt.tight_layout(0.5)
    plt.savefig('/home/marius/Dokumenter/fys4150/project3/stableearthsun-'+str(yr)+'.pgf')
    plt.savefig('/home/marius/Dokumenter/fys4150/project3/stableearthsun-'+str(yr)+'.png')



