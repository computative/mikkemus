from numpy import *
from numpy.linalg import *


class cluster:
    def __init__(self, coordinates, velocities, masses):
        self.r0 = coordinates
        self.v0 = velocities
        self.M = masses
        self.algorithm = True # true indicates verlet method
        self.n = len(coordinates)

    def method(self,algorithm):
        if (algorithm == "Euler") or (algorithm == "euler"):
            self.algorithm = False
        else:
            self.algorithm = True
        return self

    def solve(self,quantities, h, N):
        data = []
        label = []
        
        if self.algorithm:
            r, v, V = self.verlet(h,N)
        else:
            r, v, V = self.euler(h,N)
        
        for quantity in quantities:
            if quantity == 'r':
                data.append(r)
                label.append(['r'])
            elif quantity == 'v':
                data.append(v)
                label.append(['v'])
            elif quantity == 'V':
                data.append(V)
                label.append(['V'])
            elif quantity == 'T':
                print sum(v**2,axis=2)
                exit()
                data.append((0.5*sum(v**2,axis=2).T*self.M).T)
                label.append(['T'])
            elif quantity == 'p':
                data.append((r.T*self.M).T)
                label.append(['p'])
            
        return array(data), label
    
    def euler(self,h, N):
        M = self.M
        r = zeros((self.n,N,3)); v = zeros((self.n,N,3)); G = 4*pi**2;
        r[:,0] = self.r0; v[:,0] = self.v0
        F = zeros((self.n,3))

        # velocity verlet integration

        for i in arange(N-1):
            for j in arange(self.n):
                F[j] = G*sum(((r[j,i]-r[arange(0,self.n)!=j,i]).T*M[arange(0,self.n)!=j]*norm(r[j,i]-r[arange(0,self.n)!=j,i],axis=1)**-3).T,axis=0)
                r[j,i+1] = r[j,i] + h*v[j,i]
                v[j,i+1] = v[j,i] - h*F[j]
        return r,v
    
    def verlet(self,h, N):
        M = self.M
        r,v,V = zeros((self.n,N,3)), zeros((self.n,N,3)), zeros((self.n,N,3));
        G = 4*pi**2; r[:,0] = self.r0; v[:,0] = self.v0;
        F, Fpp = zeros((self.n,3)), zeros((self.n,3))

        # velocity verlet integration

        for i in arange(N-1):
            for j in arange(self.n):
                V[j,i+1] = -G*sum(M[arange(0,self.n)!=j]*norm(r[j,i]-r[arange(0,self.n)!=j,i],axis=1)**-1)
                F[j] = G*sum(((r[j,i]-r[arange(0,self.n)!=j,i]).T*M[arange(0,self.n)!=j]*norm(r[j,i]-r[arange(0,self.n)!=j,i],axis=1)**-3).T,axis=0)
                r[j,i+1] = r[j,i] + h*v[j,i] - (h**2/2.)*F[j]
            for j in arange(self.n):
                Fpp[j] = G*sum(((r[j,i+1]-r[arange(0,self.n)!=j,i+1]).T*M[arange(0,self.n)!=j]*norm(r[j,i+1]-r[arange(0,self.n)!=j,i+1],axis=1)**-3).T,axis=0)
                v[j,i+1] = v[j,i] - (h/2.)*(Fpp[j] + F[j])
        return r,v,V

# data from NASA
coords= array([
[2.915436198756346E-03, -7.544693723556148E-04, -1.390618528070157E-04],
[1.302847368172481E-01,  2.801140182961619E-01,  1.112470944910341E-02],
[7.159985056451860E-01, -1.365428669932002E-01, -4.315142412583890E-02],
[-4.798395768079091E-01,  8.566717818936551E-01, -1.684540588602532E-04],
[1.394306114932561E+00,  1.120562027349614E-02, -3.403925187984681E-02],
[-3.829376928238765E+00,  3.698200560297282E+00,  7.025249031036523E-02],
[-5.321831187579115E+00, -8.407539750328567E+00,  3.579867486469124E-01],
[1.928628493224007E+01,  5.319099485612918E+00, -2.301047163875617E-01],
[2.755280668861852E+01, -1.178859447537265E+01, -3.922183209157875E-01],
[7.459876081527229E+00, -3.191852465184250E+01,  1.257645952904755E+00]
])

velocities = 365*array([
[3.857139744374045E-06,  5.364068584843643E-06, -9.634053032378871E-08],
[-3.124383002988246E-02,  1.271888259134524E-02,  3.905552464808973E-03],
[3.688101111931239E-03,  1.978416185656664E-02,  5.837155575794271E-05],
[-1.527454798742646E-02, -8.503390397281423E-03,  3.521735189775110E-07],
[4.219149222569453E-04,  1.519504133934046E-02,  3.079148499996481E-04],
[-5.332206086068748E-03, -5.072399242009707E-03,  1.403903802689404E-04],
[4.409642240048768E-03, -2.999859826994221E-03, -1.231805785411068E-04],
[-1.074565584179083E-03,  3.608247964124384E-03,  2.733105568624108E-05],
[1.214183902132763E-03,  2.905120905745809E-03, -8.758655862606478E-05],
[3.121097822894519E-03,  8.315704139744164E-05, -8.994834256196322E-04]
])
# sun, mercury,venus,earth,mars, jupiter, saturn, uranus, nepute, pluto
masses = asarray([1, 3.285e23/1.98855e30, 4.867E24/1.98855e30, 3.003e-6, 6.39E23/1.98855e30, 1047.**-1, 5.683E26/1.98855e30, 8.681E25/1.98855e30, 1.024E26/1.98855e30, 1.31e22/1.98855e30])


# initiate solar system
h = 1e-2
N = int(300/h)
solarsystem = cluster(coords, velocities, masses)

# find the velocity and position of all the planets
data,labels = solarsystem.method("verlet").solve(['v','r'], h, N)



import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



# plot outer solar system
f = plt.figure(figsize=(0.8*4,0.8*3))
ax = plt.gca(projection='3d')
for r, j in zip(data[1],range(1,len(data[1])+1)):
    ln1=ax.plot(r[::100,0], r[::100,1], r[::100,2],'k')
    if j in [1,6,7,8,9,10]:
        ax.scatter(r[N-1:,0], r[N-1:,1], r[N-1:,2],c='k')

plt.locator_params(nbins=4)
title = r'Orbitals of the solar system to scale'
ax.set_title(title, fontsize=11, loc='center')
ax.set_xlabel(ur'$x$   [AU]', fontsize=11)
ax.set_ylabel(ur'$y$   [AU]', fontsize=11)
ax.set_zlabel(ur'$z$   [AU]', fontsize=11)
plt.tight_layout(0.5)
plt.savefig('../benchmark/outer.png')



# plot inner solar system
f = plt.figure(figsize=(0.8*4,0.8*3))
ax = plt.gca(projection='3d')
for r,j in zip(data[1,0:6],range(1,6)):
    if j == 2:
        ax.plot(r[:50:,0], r[:50,1], r[:50,2],'k')
        ax.scatter(r[50,0], r[50,1], r[50,2],c='k')
    else:
        ax.plot(r[:1000:,0], r[:1000,1], r[:1000,2],'k')
        ax.scatter(r[1000,0], r[1000,1], r[1000,2],c='k')

plt.locator_params(nbins=3)
title = r'Orbitals of inner solar system'
ax.set_title(title, fontsize=11, loc='center')
ax.set_xlabel(ur'$x$   [AU]', fontsize=11)
ax.set_ylabel(ur'$y$   [AU]', fontsize=11)
ax.set_zlabel(ur'$z$   [AU]', fontsize=11)
ax.set_zlim(-1,1)
plt.tight_layout(0.5)
plt.savefig('../benchmark/inner.png')
