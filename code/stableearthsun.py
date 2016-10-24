# -*- coding: utf-8 -*-
from numpy import *
from numpy.linalg import *
import matplotlib.pyplot as plt

# for each year in [1, 1000] we produce figures of the error
for yr,j in zip([1,100],[-3,-2]):
    
    f, ax = plt.subplots(figsize=(0.8*4,0.8*3))
    H = logspace(j,0,40)
    data = []; data2 = [] # data matrices store data
    
    # we evaluate the error for each steplength h in H
    for h in H:
        N = int(ceil(yr/h))
        
        # initiate euler-solver
        r = zeros((N,2)); v = zeros((N,2)); G = 4*pi**2;
        M = 1; r[0] = [1,0]; v[0] = (G*M)**0.5*array([-r[0,1], r[0,0]])
        for i in arange(N-1):
            r[i+1] = r[i] + h*v[i]
            v[i+1] = v[i] - h*G*M*r[i]/norm(r[i])**3
    
        # initiate verlet-solver
        r2 = zeros((N,2)); v2 = zeros((N,2));r2[0] = [1,0];
        v2[0] = (G*M)**0.5*array([-r[0,1], r[0,0]])
        for i in range(N-1):
            r2[i+1] = r2[i] + h*v2[i] - (h**2/2)*G*M*r2[i]/norm(r2[i])**3
            v2[i+1] = v2[i] - G*M*(h/2)*(r2[i+1]/norm(r2[i+1])**3 \
                                         + r2[i]/norm(r2[i])**3 )

        # we store the data for each year, steplength and method
        data.append( max( norm(r, axis=1)-1 ) )
        data2.append( max( norm(r2, axis=1)-1 ) )

    # we plot the error as a function of steplength for both methods
    ax.loglog(H,data,'k', label='Forward Euler')
    ax.loglog(H,data2,'k--', label='Velocity Verlet')
    ax.set_xlabel(ur'$\Delta t$    [yrs]')
    ax.set_ylabel(ur'Error $\|\vec{r} - \vec{s}\|_\infty$    [AU]')
    ax.set_title('Error bound after '+str(yr)+' yrs simulation', fontsize=11)
    plt.tight_layout(0.5)
    plt.savefig('../benchmark/stableearthsun-'\
                +str(yr)+'.png')
