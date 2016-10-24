# -*- coding: utf-8 -*-
from numpy import *
from numpy.linalg import *
import matplotlib.pyplot as plt


# constants and initial conditions
h = 1e-4; N = int(ceil(1/h))
r = zeros((N,2)); v = zeros((N,2)); U = zeros(N); G = 4*pi**2;
M = 1; r[0] = [1,0]; v[0] = (G*M)**0.5*array([-r[0,1], r[0,0]])
U[0] = -G*1*3.003e-6/norm(r[0])


# forward euler integration
for i in arange(N-1):
    r[i+1] = r[i] + h*v[i]
    U[i+1] =-G*1*3.003e-6/norm(r[i+1]) # storing potential energy
    v[i+1] = v[i] - h*G*M*r[i]/norm(r[i])**3

    
# velocity verlet integration
r2 = zeros((N,2)); v2 = zeros((N,2)); r2[0] = [1,0];
v2[0] = (G*M)**0.5*array([-r2[0,1], r2[0,0]]); U2 = zeros(N); 
U2[0] = -G*1*3.003e-6/norm(r2[0])
for i in range(N-1):
    r2[i+1] = r2[i] + h*v2[i] - (h**2/2)*G*M*r2[i]/norm(r2[i])**3
    U2[i+1] =-G*1*3.003e-6/norm(r2[i+1]) # storing potential energy
    v2[i+1] = v2[i] - G*M*(h/2)*(r2[i+1]/norm(r2[i+1])**3 + r2[i]/norm(r2[i]) )


# computing potential energy an momentum for both algorithms
T = 0.5*3.003e-6*norm(v,axis=1)**2
T2 = 0.5*3.003e-6*norm(v2,axis=1)**2
p = 3.003e-6*norm(v,axis=1)
p2 = 3.003e-6*norm(v2,axis=1)
t = linspace(0,1,N)

# printing the scaling factors for the figure caption
print 'T0,U0,p0: ',T[0],U[0],p[0] 

# plotting quantities for the euler-solver
f, ax = plt.subplots(figsize=(0.8*4,0.8*3))
ax.plot(t,T/T[0]-1,'k', label=ur"$T/T_0 - 1$")
ax.plot(t,U/U[0]-1,'k--', label=ur"$U/U_0 - 1$",linewidth=1.3)
ax.plot(t,p/p[0]-1,'k-.', label=ur"$p/p_0 - 1$", linewidth=2)
ax.set_title(ur'Variation $p,U$ and $T$', fontsize=11);
ax.set_xlabel(ur'$t$    [yrs]');
ax.set_ylabel(ur'$y$-axis    [$ \ \cdot \ $]');
plt.tight_layout(0.5)
plt.savefig("../benchmark/quantearthsun-a.png")

# plotting quantities for the verlet-solver
f, ax = plt.subplots(figsize=(0.8*4,0.8*3))
ax.plot(t,T2/T2[0]-1,'k', label=ur"$T/T_0 -1$")
ax.plot(t,U2/U2[0]-1,'k--', label=ur"$U/U_0 -1$",linewidth=1.3)
ax.plot(t,p2/p2[0]-1,'k-.', label=ur"$p/p_0 -1$",linewidth=2)
ax.set_title('Variation $p,U$ and $T$', fontsize=11);
ax.set_xlabel(ur'$t$    [yrs]');
ax.set_ylabel(ur'$y$-axis    [$\ \cdot \ $]');
ax.set_ylim(-2.6e-7,0.6e-7)
plt.tight_layout(0.5)
plt.savefig("../benchmark/quantearthsun-b.png")

