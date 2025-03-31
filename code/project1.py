#! C:\Users\Vidar1\Desktop\Numdiff\MittWorkspace\.venv\Scripts\python.exe

import numpy as np
import pandas as pd
import scipy as sp
from matplotlib import pyplot as plt

def f1(A,told,uold):
    return uold*A

def lotka(A,t,u):
    a = 3
    b = 9
    c = 15
    d = 15
    dxdt = a*u[0]-b*u[0]*u[1]
    dydt = c*u[0]*u[1]-d*u[1]
    dudt = np.array([dxdt, dydt])
    return dudt

def vanDerPol(mu,t,y):
    dy1dt = y[1]
    dy2dt = mu * (1-y[0]**2)*y[1]-y[0]
    return np.array([dy1dt,dy2dt])

def H(x,y):
    a = 3
    b = 9
    c = 15
    d = 15
    return c*x+b*y-d*np.log(x)-a*np.log(y)

def RK4step(A,told, uold, h, f):
    Yprim1 = f(A,told,uold)
    Yprim2 = f(A,told+h/2,uold+h*Yprim1/2)
    Yprim3 = f(A,told+h/2,uold+h*Yprim2/2)
    Yprim4 = f(A,told+h,uold+h*Yprim3)
    yn = uold + h/6*(Yprim1+2*Yprim2+2*Yprim3+Yprim4)
    print(Yprim1,Yprim2,Yprim3,Yprim4)
    return yn

def RK4int(A, y0, t0, tf, N, f):
    tgrid = np.linspace(t0,tf,N+1)
    y_list =np.array([y0])
    yold = y0
    error_list = np.array([0])
    h = (tf-t0)/N    
    for x in range(N):
        y,e = RK34step(A,tgrid[x],yold,h,f) 
        y_list = np.append(y_list, y)
        yold = y
        err = sp.linalg.norm(y0 * np.exp(tgrid[x+1] * A) - y)
        error_list = np.append(error_list,err)
    return tgrid, y_list, error_list

def errVsh(A, y0, t0, tf, f):
    N = np.array([2**13, 2**12, 2**11, 2**10, 2**9, 2**8])
    h = (tf-t0)/N
    error = np.array([])
    for x in range(len(N)):
        t,y,err = RK4int(A, y0, t0, tf, N[x], f)
        error = np.append(error, err[-1])
    return h, error

def RK34step(A, told, uold, h, f):
    Yprim1 = f(A,told,uold)
    Yprim2 = f(A,told+h/2,uold+h*Yprim1/2)
    Yprim3 = f(A,told+h/2,uold+h*Yprim2/2)
    Yprim4 = f(A,told+h,uold+h*Yprim3)
    Zprim = f(A,told+h,uold-h*Yprim1+2*h*Yprim2)
    ynew = uold + h/6*(Yprim1+2*Yprim2+2*Yprim3+Yprim4)
    znew = uold + h/6*(Yprim1+4*Yprim2+Zprim)
    
    err = sp.linalg.norm(znew-ynew)
    return ynew, err

def newstep(tol,err, errold, hold, k):
    if errold == 0:
        errold = tol
    if err == 0:
        err = tol

    hnew = ((tol/err)**(2/(3*k))*(tol/errold)**(-1/(3*k)))*hold 
    return hnew

def adaptiveRK34(A, t0, tf, y0, tol, f):
    y_list = np.array(y0) 
    t_list = np.array([t0])
    told = t0
    yold = y0
    errold = 0
    h = np.abs(tf-t0)*tol**(1/4) / (100*(1+sp.linalg.norm(f(A,t0,y0))))
    k = 4

    while (told + h < tf):
        y,err = RK34step(A, told, yold, h, f)
        told = told + h
        y_list = np.append(y_list, y)
        t_list = np.append(t_list, told)
        h = newstep(tol,err,errold,h,k)
        yold = y
        errold = err
        
        
    h = tf-told
    y,err = RK34step(A, told, yold, h, f)
    y_list = np.append(y_list,y) 
    t_list = np.append(t_list,told+h)

    return t_list,y_list

'''
t,y,e = RK4int(-2,1,0,10,2**5,f1)
plt.plot(t,y)
plt.ylabel('Y approximation')
plt.xlabel('Time')
plt.show()
h,error = errVsh(-2,1,0,10,f1)
print(h)
print(error)
plt.loglog(h,10**8*error)
plt.loglog(h,h**4)
plt.ylabel('global error')
plt.xlabel('stepsize')
plt.show()
'''
'''
y1 = RK4step(-1,0,1,0.5,f1)
y2,err = RK34step(-1,0,1,0.5,f1)
print(y1)
print(y2)
'''
'''
t,y = adaptiveRK34(-1,0,10,1,10**-6, f1)
plt.plot(t,y)
plt.show()
print(t)
print(y)
'''
'''
#print(RK4step(-1,0,1,0.5,f1))
#y = plt.plot(y)
#plt.show()
'''
'''
t,u = adaptiveRK34(-1,0,500,np.array([1/2,1/3]),10**(-6), lotka)
t1,u2 = adaptiveRK34(-1,0,500,np.array([1,1/2]),10**(-6), lotka)
x = u[range(0,len(u),2)]
y = u[range(1,len(u),2)]
x1 = u2[range(0,len(u2),2)]
y1 = u2[range(1,len(u2),2)]

#Plotta population vs time
plt.figure
plt.title("Population vs time")
plt.plot(t, x)
plt.plot(t, y)
plt.xlabel('Time')
plt.ylabel('Population')
plt.legend(["Prey","Predator"])
plt.show()

#Plotta phase function
plt.figure
plt.plot(x,y)
plt.plot(x1,y1)
plt.xlabel("Prey")
plt.ylabel("Predator")
plt.show()
'''

'''
#Test med H
t,u = adaptiveRK34(-1,0,500,np.array([1,1]),10**(-6), lotka)
x = u[range(0,len(u),2)]
y = u[range(1,len(u),2)]
H = np.absolute(H(x,y)/H(x[0],y[0])-1)
plt.plot(t,H)
plt.semilogx(t,H)
plt.title("H(x,y) as a function of time")
plt.ylabel("H(x,y)/H(x0,y0)-1")
plt.show()
'''
'''
t,y = adaptiveRK34(-1, 0, 200, np.array([2,0]),10**(-6),vanDerPol)
y1 = y[range(0,len(y),2)]
y2 = y[range(1,len(y),2)]
plt.plot(t,y2)
#plt.plot(y1,y2)
plt.title("y2 over time")
plt.xlabel("time")
plt.ylabel("y2")
plt.show()
'''
#u = np.array([1/3,1/2])
#print(x[0])
#print(x[1])
#print(u)
#print(u.item(1))
#print(lotka(1,0,u))

'''
mu = np.array([10, 15, 22, 33, 47, 68, 100, 150, 220, 330, 470, 680])
tf = 0.7*mu
steps = np.array([])
for x in range(len(mu)):
    t,y = adaptiveRK34(mu[x], 0, tf[x], np.array([2, 0]), 10**(-8), vanDerPol)
    steps = np.append(steps,len(t))

plt.loglog(mu,steps)
plt.show()
print((steps[-1]-steps[0])/(mu[-1]-mu[0]))
'''

'''
mu = np.array([10, 15, 22, 33, 47, 68, 100, 150, 220, 330, 470, 680, 1000, 10000])
tf = 0.7*mu 
t0 = 0
y = np.array([2,0])
steps = np.array([])


for x in range(len(mu)):
   t_span = np.array([t0, tf[x]])

   def vanDerPol2(t,y):
    dy1dt = y[1]
    dy2dt = mu[x] * (1-y[0]**2)*y[1]-y[0]
    return np.array([dy1dt,dy2dt])

   sol = sp.integrate.solve_ivp(vanDerPol2, y, t_span, method = 'BDF')
   steps = np.append(steps,len(sol.t))
plt.loglog(mu,steps)
plt.xlabel('mu')
plt.ylabel('Number of steps')
plt.show()
'''