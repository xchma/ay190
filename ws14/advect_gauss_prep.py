import sys,math
import numpy as np
import matplotlib.pyplot as mpl
import scipy as sp
from matplotlib import rc

L = 100

# construct the function for input initial conditions
def ini(x):
    if ((x>=0) & (x<=L)):
        return (1.0/8.0)*np.sin(2*np.pi*x/L)
    else: return 0

initial = np.vectorize(ini)

# set up the grid here. Use a decent number of zones;
# perhaps to get a dx of 0.1
x = np.arange(0,100,1)
# parameters
dx = x[1]-x[0]
v = 0.1

n = len(x)
y = np.zeros(n)
cfl = 1.0
dt = 0.3
t = 0.0

# for initial data
sigma = np.sqrt(15.0)
x0 = 30.0

#set up initial conditions
y = initial(x)

# evolve (and show evolution)
mpl.ion()
mpl.figure()
mpl.plot(x,y,'x-') # numerical data
mpl.show()

yold2 = y
yold = y
ntmax = 500

for it in range(ntmax):
    t = t+dt
    # save previous and previous previous data
    yold2 = yold
    yold = y

    method = 'upwind'
    
    # get new data; ideally just call a function
    #y = ????
    if (method=='FTCS'):
        y[1:-1] = yold[1:-1] - (0.5*yold[1:-1]*dt/dx)*(yold[2:]-yold[0:-2])
        y[0]=y[1]
        y[-1]=y[-2]
    elif (method=='upwind'):
        y[1:] = yold[1:] - (yold[1:]*dt/dx)*(yold[1:]-yold[0:-1])
        y[0]=y[1]
    elif (method=='LaxF'):
        y[1:-1] = 0.5*(yold[2:]+yold[0:-2]) \
            - (0.5*yold[1:-1]*dt/dx)*(yold[2:]-yold[0:-2])
        y[0]=y[1]
        y[-1]=y[-2]
    elif (method=='frog'):
        y[1:-1] = yold2[1:-1] - (0.5*yold[1:-1]*dt/dx)*(yold[2:]-yold[0:-2])
        y[0]=y[1]
        y[-1]=y[-2]
    elif (method=='LaxW'):
        y[1:-1] = yold[1:-1] - (0.5*yold[1:-1]*dt/dx)*(yold[2:]-yold[0:-2]) \
            + 0.5*(v**2)*(dt**2)/(dx**2)*(yold[0:-2]+yold[2:]-2*yold[1:-1])
        y[0]=y[1]
        y[-1]=y[-2]
    else:
        print "Must provide a method"
        break

    # print "it = ",it,err
    mpl.clf()
    # plot numerical result
    mpl.plot(x,y,'b-')
    mpl.xlabel('x')
    mpl.ylabel('u')
    mpl.text(1,-0.14,"t="+str(t))
    mpl.ylim([-0.15,0.15])
    mpl.draw()

mpl.show()

