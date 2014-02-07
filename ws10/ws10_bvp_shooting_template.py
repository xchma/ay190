#!/usr/bin/env python

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

def test(npoints):
    # set up grid
    xmin = 0.0
    xmax = 1.0
    # set up grid
    x = np.linspace(xmin,xmax,npoints)
    # dx based on x[1] and x[0]
    dx = x[1]-x[0]

    # boundary values
    A = 0 # inner boundary
    B = 0.1 # outer boundary

    def calc_rhs(u,xx):
        # rhs routine
        # rhs[0] is rhs for y
        # rhs[1] is rhs for u
        rhs = np.zeros(2)
        rhs[0] = u
        rhs[1] = 12*xx - 4

        return rhs

    def integrate_FE(z,x):
        # forward-Euler integrator
    
        # make an array for all points
        # entry 0 contains y
        # entry 1 contains y'
        yy = np.zeros((npoints,2))

        yy[0,0] = A # boundary value A for y at x=0
        yy[0,1] = z # guessed boundary value for y' at x=0

        for i in range(npoints-1):
            yy[i+1,:] = yy[i,:] + dx*calc_rhs(yy[i,1],x[i])

        return yy

    # get initial guess for derivative
    z0 = -1100000.0
    z1 = 10000000.0
    yy0 = integrate_FE(z0,x)
    yy1 = integrate_FE(z1,x)
    phi0 = yy0[npoints-1,0] - B
    phi1 = yy1[npoints-1,0] - B
    dphidz = (phi1-phi0)/(z1-z0) # dphi/dz

    i = 0
    itmax = 100
    err = 1.0e99
    criterion = 1.0e-12

    z0 = z1
    phi0 = phi1
    while (err > 1.0e-12 and i < itmax):
        z1 = z0 - phi0/dphidz # secand update
        yy = integrate_FE(z1,x)
        phi1 = yy[npoints-1,0] - B
        dphidz = (phi1-phi0)/(z1-z0) # dphi/dz numerical
        err = np.abs(z1-z0) # your error measure
        z0 = z1
        phi0 = phi1
        i = i+1

        print i,z1,phi1

    return x,yy[:,0]

x1,y1=test(11)
x2,y2=test(101)
x3,y3=test(1001)

xinf=np.linspace(0,1,10001)

plt.plot(x1,y1,'b-',x2,y2,'g-',x3,y3,'c-')
plt.plot(xinf,2.0*xinf**3-2*xinf**2+0.1*xinf,"r-")
plt.xlabel('x')
plt.ylabel('y')
plt.title('Shooting method')
plt.legend(["10 points","100 points","1000 points","true solution"],loc='best')
plt.show()




