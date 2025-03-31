import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

#1.1. retrieve code from P2 in order to create Tdx matrix
#and eulerstep from P1.

def G(N, L):
    xgrid = np.linspace(0, L, N)
    #Gx = [(np.sin(np.pi*x))**2 for x in xgrid]
    Gx = [np.exp(-100*(x-0.5)**2) for x in xgrid]
    #Gx = [(np.cos(np.pi*x) -1 + 2*x) for x in xgrid]
    #G = np.diag(Gx)
    print('Gx ' + str(len(Gx)))
    return Gx

def T(L, N):
    T = np.zeros((N, N))
    np.fill_diagonal(T, -2)
    np.fill_diagonal(T[:, 1:], 1)
    np.fill_diagonal(T[1:,:], 1)

    #periodic boundary conditions
    T[N-1,0]=1
    T[0,N-1]=1    

    return T

def S(L, N):
    S = np.zeros((N, N))
    #np.fill_diagonal(S, 1)
    np.fill_diagonal(S[:, 1:], 1)
    np.fill_diagonal(S[1:,:], -1)

    #periodic boundary conditions
    S[N-1,0]=1
    S[0,N-1]=-1

    return S

def eulerstep(Tdx, uold, dt):
    unew = uold + dt * np.dot(Tdx, uold)
    return unew

def solve(Tdx, gx, dt, M):
    #gx: initial cond.
    #dt: step size
    #M: number of time steps

    u = np.zeros((M+1, len(gx)))
    #initial condition
    u[0] = gx

    for i in range(1,M+1):
        u[i] = eulerstep(Tdx, u[i-1], dt)
        #u[i] = TRstep(Tdx, u[i-1],dt)

    #adding u(t,0) = u(t,L) = 0
    u0 = np.zeros((u.shape[0], 1))
    u = np.concatenate((u, u0), axis = 1)
    u = np.concatenate((u0,u), axis = 1)
    #print(u)
    return u


#1.2
def TRstep(Tdx, uold, dt):

    A = np.eye(len(uold)) - 0.5 * dt * Tdx
    b = np.dot((np.eye(len(uold)) + 0.5 * dt * Tdx), uold)
    unew = sp.linalg.solve(A,b)

    return unew


#2.1
#takes a step of size dt starting from Initial cond. u
#amu = a*dt / dx
def LaxWen(u, amu, N):
    L = np.zeros((len(u), len(u)))
    np.fill_diagonal(L, 1-amu**2)
    np.fill_diagonal(L[:, 1:], amu*0.5 * (1+amu)) #superdiagonal
    np.fill_diagonal(L[1:,:], -amu*0.5*(1-amu)) #subdiagonal
    L[0][-1] = -amu*0.5*(1-amu)
    L[-1][0] = amu*0.5 * (1+amu)
    unew = np.dot(L,u)
    RMS = sp.linalg.norm(unew)/np.sqrt(N)

    return unew, RMS


def solveLax(amu, gx, M, N):
    u = np.zeros((M+1, len(gx)))
    #initial condition
    u[0] = gx
    RMS = np.zeros(M+1)
    RMS[0] = sp.linalg.norm(gx)/np.sqrt(N)

    for i in range(1,M+1):
        u[i], RMS[i] = LaxWen(u[i-1], amu,N)

    #adding u(t,0) = u(t,L) = 0
    u0 = np.zeros((u.shape[0], 1))
    for i in range(u.shape[0]):
        u0[i] = u[i][0]

    u = np.concatenate((u, u0), axis = 1)
    #u = np.concatenate((u0,u), axis = 1)
    print('u rows: ' + str(len(u)) + ', cols: ' + str(len(u[0])))
    

    return u, RMS


#3.1
def convdif(u, a, d, dt, dx, L, N):
    #a = convection velocity
    #d = diffusivity
    #u function

    Sdx = S(L,N) / (2*dx)
    Tdx = T(L,N) / (dx**2)

    A = d*Tdx - a*Sdx
    AA = np.eye(N,N) - dt*0.5*A
    B = np.dot(u, np.eye(N,N) + dt*0.5*A)
    
    # u2 = sp.linalg.solve(AA,B)
    u2 = np.dot(B, np.linalg.inv(AA))
    print(u2.shape)
    return u2


def solveConv(gx, a, d, dt, dx, L, N, M):
    #gx: initial cond.
    #dt: step size
    #M: number of time steps

    u = np.zeros((M+1, len(gx)))
    #initial condition
    u[0] = gx
    for i in range(1,M+1):
        u[i] = convdif(u[i-1], a, d, dt, dx, L, N)

    u0 = np.zeros((u.shape[0], 1))

    for i in range(u.shape[0]):
        u0[i] = u[i][0]

    u = np.concatenate((u, u0), axis = 1)

    return u