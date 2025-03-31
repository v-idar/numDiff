import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from functions3 import G, T, eulerstep, solve, TRstep, LaxSolve

'''1.1
L, tEnd, N, M = 1, 1, 100, 200 # Längd, sluttid, Antal x-punkter och antal t-punkter
dx, dt = L/(N+1), tEnd/M
gx = G(N, L)
CFL = dt / (dx**2)
print(CFL)
Tdx = T(L, N) / (dx**2)
xx = np.linspace(0, L, N+2)
tt = np.linspace(0, tEnd, M+1)

T, X = np.meshgrid(tt,xx)
print('T ' + str(len(T)))
print('X ' + str(len(X)))


approx = solve(Tdx, gx, dt, M, TRstep)
approx = approx.T
print('approx, rows x cols: ' + str(len(approx)) + ' x ' + str(len(approx[0])))

# Plotting the 3D surface
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(T, X, approx, cmap='viridis')
ax.set_xlabel('Time (t)')
ax.set_ylabel('Position (x)')
ax.set_zlabel('Solution (u)')
ax.set_title('Numerical Solution u(t, x)')
plt.colorbar(surf, shrink=0.5, aspect=5)  # Colorbar
plt.show()
'''

#2.1
#'''
a = 1.5
L = 1
tEnd = 5
N = 100
M = 500
dx = L/N
dt = tEnd/M
print('dx ' + str(dx))
print('dt ' + str(dt))

amu = a*dt/dx
gx = G(N,L)
CFL = amu

xx = np.linspace(0,L,N+1)
tt = np.linspace(0,tEnd, M+1)
T, X = np.meshgrid(tt,xx)
print('T ' + str(len(T)))
print('X ' + str(len(X)))
print('CFL : ' + str(CFL))
lax, RMS = LaxSolve(amu, gx, M,N)
lax = lax.T
# '''


# Plotting the 3D surface
'''
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(T, X, lax, cmap='viridis')
ax.set_xlabel('Time (t)')
ax.set_ylabel('Position (x)')
ax.set_zlabel('Solution (u)')
ax.set_title('Numerical Solution u(t, x)')
plt.colorbar(surf, shrink=0.5, aspect=5)  # Colorbar
plt.show()

plt.plot(tt, RMS)
plt.show()
'''


#3.1
'''
N = 50
M = 500
L = 1
tEnd = 1
a = 1
d = 100
dx = L/N
dt = tEnd/M
Pe = np.abs(a/d)
gx = G(N,L)
print('Pe: ' + str(Pe))
approx = solveConv(gx,a,d,dt, dx, L, N, M)
approx = approx.T
print('U: rows: ' +  str(len(approx)) + ', cols: ' + str(len(approx[0])))

xx = np.linspace(0, L, N+1)
tt = np.linspace(0, tEnd, M+1)
T, X = np.meshgrid(tt, xx)
print('T ' + str(len(T)))
print('X ' + str(len(X)))


fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(T, X, approx, cmap='viridis')
ax.set_xlabel('Time (t)')
ax.set_ylabel('Position (x)')
ax.set_zlabel('Solution (u)')
ax.set_title('Numerical Solution u(t, x)')
plt.colorbar(surf, shrink=0.5, aspect=5)  # Colorbar
plt.show()
'''