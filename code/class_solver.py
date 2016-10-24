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
                # endre kreftene dersom sola er stasjonaer
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