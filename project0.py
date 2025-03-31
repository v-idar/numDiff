import numpy as np
import pandas as pd
import scipy as sp
from matplotlib import pyplot as plt

#Task 1.1
def unew(A,uold, h):
    return uold + h*A*uold

#Task 1.2/1.4
def eulerint(A, y0, t0, tf, N):
    
    y = y0
    y_list =np.array([y0])
    error_list = np.array([0])
    tgrid = np.linspace(t0,tf,N+1)
    h = (tf-t0)/N
    
    for x in range(N):
        y = unew(A,y,h)
        y_list = np.append(y_list, y)
        err = np.linalg.norm(y0*sp.linalg.expm(tgrid[x]*A)-y)
        error_list = np.append(error_list,err)
    return tgrid, y_list, error_list


#Task 1.3/1.4
def errVsh(A, y0, t0, tf):
    N = np.linspace(1,1000,10,dtype=int)
    global_error = np.array([])
    for x in range(len(N)):
        t,a_list,e_list = eulerint(A,y0,t0,tf,N[x])
        global_error = np.append(global_error,e_list[-1])
    return N,global_error



#Test1.3
'''
N, error_list = errVsh(5,1,1,2)
error_list_norm = np.linalg.norm(error_list)
plt.plot(N,error_list)
plt.ylabel('global error')
plt.xlabel('number of steps')
plt.yscale('log')
plt.xscale('log')
plt.show()
'''
#Test1.4
'''
time, approx, error = eulerint(-1,1,0,2,1000)
plt.figure()
plt.plot(time,error)
plt.xlabel("Time")
plt.ylabel("Error")
plt.yscale('log')
plt.show()
'''
#Test relative error
'''
print(approx.shape)
rel_error = error/approx
plt.figure()
plt.plot(time,rel_error)
plt.xlabel("Time")
plt.ylabel("Relative Error")
plt.yscale('log')
plt.show()
'''
#Test final error att different stepsizes
x = 2**np.arange(-5,-10,-1)
#for i in range(5):