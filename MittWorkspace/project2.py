import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from sklearn import preprocessing
from functions import errVsh, f, twopBVP, eigen, twopBVP_SL, errVsNSL, eigenWave, V

# # Task 1.1
# alpha, beta, L, N = 0, 0, np.pi, 100 # BC, length and Steps
# fvec = f(alpha, beta, L, N)         # Second derivatives
# y, yreal, tgrid, err = twopBVP(fvec, alpha, beta, L, N)
# dx, errorlist = errVsh(alpha, beta, L)

# # Plot expexted twopBVP and expected y solution over time
# plt.plot(tgrid, y, label='twopBVP solution')
# plt.plot(tgrid, yreal, label= 'expected solution')
# plt.xlabel('time')
# plt.ylabel('solution, Y')
# plt.title('two point BVP')
# plt.legend()
# plt.grid(True)
# plt.show()


'''# Plot error over stepsize
plt.loglog(dx, errorlist, label='error')
plt.loglog(dx, dx**2, label = 'expected error') 
plt.xlabel('step size h')
plt.ylabel('error')
plt.title('error diagram')
plt.legend()
plt.grid(True)
plt.show()
'''

'''# Task 1.2
E, q = 1.9*10**11, -50*10**3 # Young's modulus and load density
alpha, beta, L, N = 0, 0, 10, 999 # Start parameters
dx = L/(N+1)    # Stepsize
xgrid = np.arange(dx, L, dx) # Grid
I = 10**(-3)*(3-2*(np.cos(np.pi*xgrid/L))**12) # Crossection moment of inertia
Mbiss = q * np.ones(N)
M, Mreal, x, errM = twopBVP(Mbiss, alpha, beta, L, N)
ubiss = M[1:-1] / (E*I)
u, ureal, x, erru = twopBVP(ubiss, alpha, beta, L, N)
print(u[500])
'''

'''# Plot the computed deflection
plt.plot(x, u, label='twopBVP solution')
plt.xlabel('x')
plt.ylabel('solution, u')
plt.title('The beam equation')
plt.legend()
plt.grid(True)
plt.show()
'''

# Task 2.1 
alpha, beta, L, N = 0, 0, 1, 499

# Task 2.1
Nvec = [2**5, 2**6, 2**7, 2**8, 2**9, 2**10]
for i in range(3):
    eigreal = -np.pi**2/4 * (2*i+1)**2
    errlist = errVsNSL(alpha, beta, L, Nvec, eigreal, i)
    plt.loglog(Nvec, errlist, label='error' + str(i))
plt.loglog(Nvec, [n**(-2) for n in Nvec], label= 'expected error')
plt.xlabel('N')
plt.ylabel('error')
plt.title('errVsh')
plt.legend()
plt.grid(True)
plt.show()


# eigvalues, eigvectors = eigen(alpha, beta, L, N)    
# dx = L/N
# xgrid = np.linspace(0, L, N+1)
# for i in range(3):
#     eigvector = np.insert(eigvectors[:,i], 0, 0)
#     plt.plot(xgrid, eigvector, label= 'eigenfunction ' + str(i+1))
# plt.xlabel('grid')
# plt.ylabel('eigenfunctions')
# plt.title('Eigen function')
# plt.legend()
# plt.grid(True)
# plt.show()




# Test 2.2 Plottar dem första vågfunktionerna samt energinivåer

# dx = L / (N+1)
# v = 0
# xgrid = np.linspace(0, L, N+1)
# print(xgrid)
# n = 3
# eigvalues, eigvectors = eigenWave(alpha, beta, L, N, n, 0)


# # # Waveequations
# # for i in range(n):
# #     eigvector = np.insert(eigvectors[:,i], 0, 0)
# #     eigvector = np.append(eigvector, 0)
# #     E = v-eigvalues[N-2-i]
# #     eigvector = eigvector / np.linalg.norm(eigvector)   
# #     plt.axhline(y=E)
# #     plt.plot(xgrid, eigvector + E, label= 'eigenfunction ' + str(i+1))
# # plt.show()

# # Probdensities
# for i in range(n):
#     eigvector = np.insert(eigvectors[:,i], 0, 0)
#     eigvector = np.append(eigvector, 0)
#     probdensities = [eig**2 for eig in eigvector]
#     E = v-eigvalues[N-2-i]
#     probdensities = probdensities / np.linalg.norm(probdensities)   
#     plt.axhline(y=E)
#     plt.plot(xgrid, probdensities + E, label= 'Probability densities' + str(i+1))
# plt.show()

# Task 2.2: non zero potentials
# N = 499
# L = 1
# n = 6
# alpha = beta = 0

# P = V(N,L)
# xgrid = np.linspace(0, L, N+1)
# eigvalues, eigvectors = eigenWave(alpha, beta, L, N, n, P)


# # Waveequations
# for i in range(n):
#      eigvector = np.insert(eigvectors[:,i], 0, 0)
#      eigvector = np.append(eigvector, 0)
#      eigvector = eigvector*600
#      E = P[i][i]-eigvalues[N-2-i]
#      #eigvector = eigvector / np.linalg.norm(eigvector) 
#      #probdensities = [eig**2 for eig in eigvector]
#      plt.axhline(y=E)
#      plt.plot(xgrid, eigvector + E, label= 'eigenfunction ' + str(i+1))
# plt.xlabel('X')
# plt.ylabel('eigen functions')
# plt.title('tunneling in potential well')
# plt.legend()
# plt.show()