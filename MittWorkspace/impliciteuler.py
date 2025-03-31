#! C:\Users\Vidar1\Desktop\Numdiff\MittWorkspace\.venv\Scripts\python.exe

import numpy as np
import pandas as pd
import scipy as sp
from matplotlib import pyplot as plt

#Task 1.1
def ieulerstep(A,uold, h):
    return uold*(1-h*A)**-1

#Task 1.2/1.4
def ieulerint(A, y0, t0, tf, N):
    
    y = y0
    y_list =np.array([y0])
    error_list = np.array([0])
    tgrid = np.linspace(t0,tf,N+1)
    h = (tf-t0)/N
    
    for x in range(N):
        y = ieulerstep(A,y,h)
        y_list = np.append(y_list, y)
        err = sp.linalg.norm(y0*sp.linalg.expm(tgrid[x]*A)-y)
        error_list = np.append(error_list,err)

    return tgrid, y_list, error_list


#Task 1.3/1.4
def ierrVsh(A, y0, t0, tf):
    N = 2**np.arange(5,15,1,dtype=int)
    global_error = np.array([])
    for x in range(len(N)):
        t,a_list,e_list = ieulerint(A,y0,t0,tf,N[x])
        global_error = np.append(global_error,e_list[-1])
    return N,global_error


'''
#Test1.2
t,y_list,error_list = eulerint(1,1,0,1,100)
error_list_norm = np.linalg.norm(error_list)
plt.plot(y_list)
plt.ylabel('approx')
plt.xlabel('time grid')
plt.yscale('log')
plt.xscale('log')
plt.show()
'''

#Test1.3
#'''
for x in np.arange(-4,5,1):
    N, error_list = ierrVsh(x,1,0,1)
    plt.plot(N,error_list,label='$\lambda$' + ' = ' + str(x))
plt.legend(fontsize = "8",loc = "lower left")
plt.ylabel('global error')
plt.xlabel('number of steps')
plt.yscale('log')
plt.xscale('log')
plt.show()
#'''
#Test1.4
'''
for x in np.arange(-4,5,1):
    time, approx, error = eulerint(x,1,0,1,2**5)
    plt.plot(time,error,label='$\lambda$' + ' = ' + str(x))
plt.legend(fontsize = "8",loc = "lower left")
plt.text(0.1,3, "N = " + str(2**5), style = 'italic', bbox={'facecolor': 'blue', 'alpha': 0.5, 'pad': 8})
plt.xlabel("Time")
plt.ylabel("Error")
plt.yscale('log')
plt.show()
'''
#Test relative error
'''
for x in np.arange(-4,5,1):
    time, approx, error = eulerint(x,1,0,100,2**10)
    rel_error = error/sp.linalg.norm(np.exp(time*x))
    plt.plot(time,rel_error,label='$\lambda$' + ' = ' + str(x))
plt.legend(fontsize = "8",loc = "lower left")
plt.text(0.1,3, "N = " + str(2**5), style = 'italic', bbox={'facecolor': 'blue', 'alpha': 0.5, 'pad': 8})
plt.xlabel("Time")
plt.ylabel("Relative error")
plt.yscale('log')
plt.show()
'''
#Test matrices
'''
x = np.array([[-1,-5],[0,-3]])
y0 = np.array([[1],[1]])
time, approx, error = eulerint(x,y0,0,10,2**10)
#rel_error = error/sp.linalg.norm(np.exp(time*x))
plt.plot(time,error,label='$\lambda$' + ' = ' + str(x))
plt.legend(fontsize = "8",loc = "lower left")
plt.text(0.1,3, "N = " + str(2**5), style = 'italic', bbox={'facecolor': 'blue', 'alpha': 0.5, 'pad': 8})
plt.xlabel("Time")
plt.ylabel("error")
plt.yscale('log')
plt.xscale('log')
plt.show()
'''
