import sys,math
import numpy as np
import matplotlib.pyplot as mpl
import scipy as sp
from matplotlib import rc

def apply_bcs(y):
    # apply boundary conditions
    # you need to fill in code
    y[0]=y[1]
    return y


def analytic(t,x):
    return np.exp(-(x-v*t-x0)**2/(2*sigma**2))

# set up the grid here. Use a decent number of zones;
# perhaps to get a dx of 0.1
x = np.arange(0,100,0.1)
# parameters
dx = x[1]-x[0]
v = 0.1

n = len(x)
y = np.zeros(n)
cfl = 1.0
dt = 0.9
t = 0.0

# for initial data
sigma = np.sqrt(15.0)
x0 = 30.0

#set up initial conditions
y = analytic(t,x)

rc('text',usetex=True)

# evolve (and show evolution)
mpl.ion()
mpl.figure()
mpl.plot(x,y,'x-') # numerical data
mpl.plot(x,analytic(t,x),'r-') # analytic data
mpl.show()

yold2 = y
yold = y
ntmax = 100

for it in range(ntmax):
    t = t+dt
    # save previous and previous previous data
    yold2 = yold
    yold = y

    method = 'FTCS'
    
    # get new data; ideally just call a function
    #y = ????
    if (method=='FTCS'):
        y[1:-1] = yold[1:-1] - (0.5*v*dt/dx)*(yold[2:]-yold[0:-2])
        y[0]=y[1]
        y[-1]=y[-2]
    elif (method=='upwind'):
        y[1:] = yold[1:] - (v*dt/dx)*(yold[1:]-yold[0:-1])
        y[0]=y[1]
    elif (method=='LaxF'):
        y[1:-1] = 0.5*(yold[2:]+yold[0:-2]) - (0.5*v*dt/dx)*(yold[2:]-yold[0:-2])
        y[0]=y[1]
        y[-1]=y[-2]
    elif (method=='frog'):
        y[1:-1] = yold2[1:-1] - (0.5*v*dt/dx)*(yold[2:]-yold[0:-2])
        y[0]=y[1]
        y[-1]=y[-2]
    elif (method=='LaxW'):
        y[1:-1] = yold[1:-1] - (0.5*v*dt/dx)*(yold[2:]-yold[0:-2]) \
            + 0.5*(v**2)*(dt**2)/(dx**2)*(yold[0:-2]+yold[2:]-2*yold[1:-1])
        y[0]=y[1]
        y[-1]=y[-2]
    else:
        print "Must provide a method"
        break


    # after update, apply boundary conditions
    # apply_bcs(x,y) 
    # y = apply_bcs(y)

    # get analytic result for time t
    yana = analytic(t,x)
    # compute error estimage
    err = 0
    # err = ???
    err = yana - y
    # print "it = ",it,err
    mpl.clf()
    # plot numerical result
    mpl.plot(x,y,'b-')
    # plot analytic results
    mpl.plot(x,yana,'r-')
    # plot error
    mpl.plot(x,err,'g-')
    mpl.ylim([-0.05,1.0])
    mpl.draw()

mpl.legend(["Numerical","Analytical","Error"],loc='best')
mpl.xlabel('x')
mpl.ylabel('u')
mpl.text(60,0.2,r"FTCS")
mpl.text(60,0.1,r"$\alpha=0.9,~t=90,~\sigma=\sqrt{15}$")
mpl.savefig("fig4.pdf")
#mpl.show()

