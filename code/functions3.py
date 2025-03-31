import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

#1.1. retrieve code from P2 in order to create Tdx matrix
#and eulerstep from P1.

def G(N, L):
    xgrid = np.linspace(0,L,N)
    Gx = [600*np.sin(np.pi*x/2) - 600*x for x in xgrid]
    #G = np.diag(Gx)
    print('Gx ' + str(len(Gx)))
    return Gx

def T(L, N):
    T = np.zeros((N, N))
    np.fill_diagonal(T, -2)
    np.fill_diagonal(T[:, 1:], 1)
    np.fill_diagonal(T[1:,:], 1)
    
    return T


def eulerstep(Tdx, uold, dt):
    # Takes an explicit eulerstep with the inparam: Matrix Tdx, old u-value
    # and timestep dt.
    unew = uold + dt * np.dot(Tdx, uold)
    return unew

def solve(Tdx, gx, dt, numT, f):
    #gx: initial cond.
    #dt: step size
    #numT: number of time steps
    
    u = np.zeros((numT+1, len(gx)))     
    u[0] = gx # Initial condition
    for i in range(1,numT+1): # Räknar ut u(i, x)
                u[i] = f(Tdx, u[i-1], dt)
    u0 = np.zeros((u.shape[0], 1)) # Boundary condition
    u = np.concatenate((u, u0), axis = 1) #adding u(t,0) = 0
    u = np.concatenate((u0,u), axis = 1) #adding u(t,L) = 0
    return u

#1.2
def TRstep(Tdx, uold, dt):
    # The Crank-Nicolson method. The same as applying the trapezoidal rule 
    # to the equation
    A = np.eye(len(uold)) - 0.5 * dt * Tdx
    b = uold + 0.5 * dt * np.dot(Tdx, uold)
    unew = sp.linalg.solve(A,b)
    return unew

def LaxWen(u, amu):
    T = np.zeros((len(u), len(u)))
    np.fill_diagonal(T, 1-amu**2)
    np.fill_diagonal(T[:, 1:], amu/2 * (1+amu))
    np.fill_diagonal(T[1:,:], -amu/2*(1-amu))

    unew = np.dot(T * u)

    return unew

def LaxSolve(amu, gx, dt, numT):
    #gx: initial cond.
    #dt: step size
    #numT: number of time steps
    u = np.zeros((numT+1, len(gx)))
    #initial condition
    u[0] = gx

    RMS = np.zeros(M+1)
    RMS[0] = sp.linalg.norm(gx)/np.sqrt(N)

    for i in range(1,M+1):
        u[i], RMS[i] = LaxWen(u[i-1], amu,N)

    for i in range(1,numT+1): # Räknar ut u(i, x)
        u[i] = LaxWen(Tdx, u[i-1], dt)
    #adding u(t,0) = u(t,L) = 0
    u0 = np.zeros((u.shape[0], 1))
    u = np.concatenate((u, u0), axis = 1)
    u = np.concatenate((u0,u), axis = 1)
    #print(u)
    return u

def convdif(u, a, d, dt, dx):
    # The single step should use the trapezoidal rule
    # Boundary conditions should  
    unew  =  0
    return unew