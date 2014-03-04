#!/usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as mpl

## This is a 1D planar SPH code

#############################################
# General constants
twothirds = 2.0/3.0

#############################################
# Problem-specific constants
gamma = 1.4

#############################################
# Function definitions

def smoothing_kernel(r,h):
    """
    Smoothing kernel
    """
    invh = 1.0/h
    u = r*invh
    alphaD = twothirds * invh
    if u >= 0.0 and u < 1.0:
        return alphaD*(1-1.5*u**2+0.75*u**3)
    elif u < 2.0:
        return alphaD*0.25*(2.0-u)**3
    else:
        return 0.0


def dsmoothing_kernel_rij(rij,h):
    """ 
    Derivative of the smoothing kernel for pairs
    d/dr_i W(rij,h)
    """

    invh = 1.0/h
    u = np.fabs(rij) * invh
    alphaDoh = twothirds * invh * invh
    if u >= 0.0 and u < 1.0:
        return alphaDoh*(3.0*rij/h*(0.75*u - 1.0))
    elif u < 2.0:
        return alphaDoh*-0.75*(2.0-u)**2 * rij/np.fabs(rij)
    else:
        return 0.0


def neighbor_search_1D(r,h,n):
    """
    For each particle i, find its neighbors that
    are closer than 2 times the smoothing length h
    """
    neibs = np.zeros( (n,n), dtype=int)
    nneibs = np.zeros(n,dtype=int)

    # inefficient n*n scaling, but straightforward
    for i in range(n):
        for j in range(n):
            if np.fabs(r[i]-r[j]) < 2.0*h:
                neibs[i,nneibs[i]] = j
                nneibs[i] += 1

    return (neibs,nneibs)

def get_density(r,h,m,n,neibs,nneibs):
    """
    Compute the smoothed density
    """
    rhoav = np.zeros(n)
    for i in range(n):
        for j in range(nneibs[i]):
            rij = np.fabs(r[neibs[i,j]]-r[i])
            rhoav[i] += m[neibs[i,j]] * smoothing_kernel(rij,h)

    return rhoav

def get_artificial_viscosity(r,v,cs,rho,h,neibs,nneibs):
    """
    Compute the artificial viscosity for each particle i
    with each of its neighbors j
    """
    av = np.zeros((n,n))
    alpha = 1.0
    beta = 1.0
    varphi = 0.1 * h

    for i in range(n):
        for j in range(nneibs[i]):
            jj = neibs[i,j]
            rij = r[i]-r[jj]
            vij = v[i]-v[jj]
            if(rij*vij < 0.0):
                phi = h * rij*vij / (rij**2 + varphi**2)
                cij = 0.5*(cs[i]+cs[jj])
                rhoij = 0.5*(rho[i]+rho[jj])
                av[i,jj] = (-alpha*cij*phi + beta*phi**2) / \
                           rhoij
                
    return av

def get_accels(r,m,p,rho,av,h,n,neibs,nneibs):
    """ 
    Compute the accelerations (RHS of momentum equation)
    """
    accels = np.zeros(n)
    for i in range(n):
        for j in range(nneibs[i]):
            jj = neibs[i,j]
            rij = r[i]-r[jj]
            dW = dsmoothing_kernel_rij(rij,h)
            accels[i] -= m[jj]*(p[i]/rho[i]+p[jj]/rho[jj]**2+av[i,jj])*dW

    return accels
                                   

def get_energyRHS(r,m,p,rho,v,av,h,n,neibs,nneibs):
    """
    Compute the RHS of the energy equation
    """
    epsrhs = np.zeros(n)
    for i in range(n):
        for j in range(nneibs[i]):
            jj = neibs[i,j]
            rij = r[i]-r[jj]
            vij = v[i]-v[jj]
            dW = dsmoothing_kernel_rij(rij,h)
            epsrhs[i] += 0.5*m[jj]*(p[i]/rho[i]+p[jj]/rho[jj]**2+av[i,jj])*vij*dW

    return epsrhs

def get_dt(h,cs,odt,cfl):
    """
    Simplemost computation of the time step
    """
    dt = cfl*h/max(cs)
    return np.fmin(1.1*odt,dt)

def set_bcs(y,nghost,value):
    """
    Set conditions at inner and outer boundary
    """
    y[0:nghost] = value
    y[-nghost:] = value

    return y

#############################################
# Main code
# keep it simple, use numpy arrays for everything

n = 1000
nghost = 10 # number of particles at the boundaries that we ignore
xmin = -0.5
xmax = 0.5

# set up positions
r = np.linspace(xmin,xmax,n)

# set up initial analytic densities
rho = np.ones(n) 
rho[r<=0.0] = 1.0
rho[r>0.0]  = 0.25

# get masses from density:
dx = r[1]-r[0]
m = dx*rho[:]

# set smoothing length
h = dx*5.0

# set specific internal energy
eps = np.zeros(n)
eps[r<=0.0] = 2.5
eps[r>0.0]  = 1.795

# get pressur from "Gamma-Law"
press = (gamma-1.0)*rho*eps
# get speed of sound at constant entropy
cs2 = (gamma-1.0)*rho*eps + press/rho**2 * (gamma-1.0)*rho
cs = np.sqrt(cs2)
# set all velocities to zero
v = np.zeros(n)
vh = np.zeros(n) # velocities at half step

# get neighbors
(neibs,nneibs) = neighbor_search_1D(r,h,n)

# get initial densities with smoothing kernel
rho = get_density(r,h,m,n,neibs,nneibs)

# initialize old timestep to something huge
odt = 10.0e0
# set CFL factor
cfl = 0.3
# get initial timestep
dt = get_dt(h,cs,odt,cfl)
# set time to zero
time = 0.0
# end after ntmax timesteps
ntmax = 500

print 0, time, dt

mpl.ion()
mpl.figure()

# time integration loop
for nt in range(1,ntmax):

    # get artificial viscosity
    av = get_artificial_viscosity(r,v,cs,rho,h,neibs,nneibs)

    # get accelerations
    accels = get_accels(r,m,press,rho,av,h,n,neibs,nneibs)
    accels = set_bcs(accels,nghost,0.0)

    # update vh from t-1/2dt to t+1/2dt
    vh = vh + dt * accels

    # get epsrhs
    epsrhs = get_energyRHS(r,m,press,rho,vh,av,h,n,neibs,nneibs)
    epsrhos = set_bcs(epsrhs,nghost,0.0)

    # update eps
    eps = eps + dt*epsrhs

    # update v to t+dt
    v = vh + 0.5*dt * accels

    # update positions
    r = r + dt * vh

    # get new neighbors
    (neibs,nneibs) = neighbor_search_1D(r,h,n)

    # update density, pressure, cs
    rho = get_density(r,h,m,n,neibs,nneibs)
    press = (gamma-1.0)*rho*eps
    cs2 = (gamma-1.0)*rho*eps + press/rho**2 * (gamma-1.0)*rho
    cs = np.sqrt(cs2)
    odt = dt
    dt = get_dt(h,cs,odt,cfl)

    # do some output
    if nt == 1 or nt % 5 == 0:
        mpl.clf()
        mpl.plot(r[nghost:len(rho)-nghost],rho[nghost:len(rho)-nghost],"+")
        mpl.xlim([-0.6,0.6])
        mpl.ylim([0.2,1.05])
        mpl.text(-0.58,0.21,"t=%f" % time)
        mpl.draw()
        
        fname = "./data/out_%04d.dat" % (nt)
        pltname = "./plot/plot_%04d.pdf" % (nt)
        outfile = open(fname,"w")
        for i in range(n):
            outstring = "%6d %15.6E %15.6E %15.6E %15.6E \n" % \
                        (i,r[i],rho[i],press[i],v[i])
            outfile.write(outstring)
        outfile.close()
        mpl.savefig(pltname,format='pdf')
            

    time = time + dt
    print nt, time, dt

mpl.ioff()    
mpl.plot(r,rho,"+")
mpl.show()


