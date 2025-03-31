import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

def twopBVP(fvec, alpha, beta, L, N):

    dx = L / (N+1) # N + 1 Steps
    tgrid = np.linspace(0, L, N+2) # N + 2 points equals N + 1 Steps
    T = np.zeros((N, N))
    np.fill_diagonal(T, -2)
    np.fill_diagonal(T[:, 1:], 1)
    np.fill_diagonal(T[1:,:], 1)
    T = T / (dx**2)

    y = np.linalg.solve(T, fvec)
    y = np.insert(y, 0, alpha) # Add startpoint BC
    y = np.append(y, beta)  # Add enpoint BC 
    yreal = np.sin(tgrid)   # Expected value
    err = sp.linalg.norm(yreal - y)*np.sqrt(N)
    return y, yreal, tgrid, err

def errVsh(alpha, beta, L):
    N = np.array([2**13, 2**12, 2**11, 2**10, 2**9, 2**8])
    dx = L/(N+1)
    error = np.array([])
    fvec = [f(alpha, beta, L, Ni) for Ni in N]
    for Ni in N:
        fvec = f(alpha, beta, L, Ni)
        y, yreal, tgrid, err = twopBVP(fvec, alpha, beta, L, Ni)
        error = np.append(error, err)
    return dx, error

def f(alpha, beta, L, N):
    tgrid = np.linspace(0, L, N+2)
    deltax = L / (N+1)
    fvec = -np.sin(tgrid[1:-1])
    fvec[0]-=alpha/(deltax**2)
    fvec[-1]-= beta/(deltax**2)
    return fvec

def twopBVP_SL(fvec, alpha, beta, L, N):
    dx = L / (N+1)
    tgrid = np.linspace(0, L, N+2)
    T = np.zeros((N+1, N+1))
    np.fill_diagonal(T, -2)
    np.fill_diagonal(T[:, 1:], 1)
    np.fill_diagonal(T[1:,:], 1)
    T[(N, N-1)] = 2
    T = T / (dx**2)

    y = np.linalg.solve(T, fvec)
    y = np.insert(y, 0, alpha)
    y = np.delete(y, -1)      # Sätt in nya yN
    err = sp.linalg.norm(np.sin(tgrid) - y)

    return y, tgrid, err


def eigen(alpha, beta, L, N):
    dx = L / N
    T = np.zeros((N, N))
    np.fill_diagonal(T, -2)
    np.fill_diagonal(T[:, 1:], 1)
    np.fill_diagonal(T[1:,:], 1)
    T[(N-1, N-2)] = 2
    T = T/(dx**2)
    eigvalues, eigvectors = sp.linalg.eig(T)
    eigvectorsort = np.zeros((N, 3))
    eigvaluessorted = np.sort(eigvalues)
      
    for i in range(3):
        idx = np.where(eigvalues==eigvaluessorted[498-i])
        eigvectorsort[:, i] = eigvectors[:, idx[0][0]]
        print(eigvalues[idx])

    return eigvalues, eigvectorsort

def errVsNSL(alpha, beta, L, Nvec, eigreal, k):
    errlist = np.array([])
    for i in range(len(Nvec)):
        eigvalues, eigvec = eigen(alpha, beta, L, Nvec[i])
        err = np.abs(eigvalues[k]-eigreal)
        print(eigvalues[k])
        errlist = np.append(errlist, err)

    return errlist

def V(N, L):
    xgrid = np.linspace(0,L,N-1)

    #Vx = [700*(0.5 - np.abs(x-0.5)) for x in xgrid]
    Vx = [800*(np.sin(np.pi*x))**2 for x in xgrid]
    
    #choice of test potential function in well
    #Vx = [700*np.exp(x) for x in xgrid]
    V = np.diag(Vx)
    return V

def eigenWave(alpha, beta, L, N, n, V):
    dx = L / N
    T = np.zeros((N-1, N-1))
    np.fill_diagonal(T, -2)
    np.fill_diagonal(T[:, 1:], 1)
    np.fill_diagonal(T[1:,:], 1)
    T = T/(dx**2) + V
    eigvalues, eigvectors = np.linalg.eig(T)
    eigvectorsort = np.zeros((N-1, n+1))
    eigvaluessorted = np.sort(eigvalues)
      
    for i in range(n):
        idx = np.where(eigvalues==eigvaluessorted[N-2-i])
        eigvectorsort[:, i] = eigvectors[:, idx[0][0]]

        print(eigvalues[idx])
    return eigvaluessorted, eigvectorsort
