# -*- coding: utf-8 -*-
from numpy import *
from numpy.linalg import *
import matplotlib.pyplot as plt


# load path of earth and jupiter
r = loadtxt("../resources/earth.txt")
s = loadtxt("../resources/jupiter.txt")
N = len(r)

# plotting result
f, ax = plt.subplots(figsize=(0.8*4,0.8*3))
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
plt.savefig("../benchmark/earthjupiter.png")
