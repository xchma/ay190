#!/usr/bin/env python

import numpy as np
import scipy as sp


# global constants
ggrav = 6.67e-8
msun  = 1.99e33

# EOS parameters
# for white dwarfs:
polyG = 4.0/3.0
polyK = 1.244e15*0.5**polyG


#######################################
# function definitions
def RHS_HE(rad,p,rho,m):
    return -ggrav*m*rho/rad**2

def RHS_MC(rad,p,rho,m):
    return 4*np.pi*rho*rad**2

def press_to_rho(press):
    return (press/polyK)**(1.0/polyG)

def tov_RHS(rad,p,rho,m):
    
    # RHS function
    
    rhs = np.zeros(2)
    if(rad > 1.0e-10):
        rhs[0] = RHS_HE(rad,p,rho,m)
        rhs[1] = RHS_MC(rad,p,rho,m)
    else:
        rhs[0] = 0.0
        rhs[1] = 0.0

    return rhs

def tov_integrate_FE(rad,dr,p,rho,m):

    # Forward-Euler Integrator

    new = np.zeros(2)
    old = np.zeros(2)
    old[0] = p
    old[1] = m

    # forward Euler integrator
    new = old + tov_RHS(rad,p,rho,m)*dr
    
    # assign outputs
    pnew = new[0]
    mnew = new[1]

    return (pnew,mnew)

def tov_integrate_RK2(rad,dr,p,rho,m):
    
    k1 = np.zeros(2)

    k1[0] = RHS_HE(rad,p,rho,m)*dr
    k1[1] = RHS_MC(rad,p,rho,m)*dr

    rho_plus_half_k1 = press_to_rho(p+0.5*k1[0])

    k2 = np.zeros(2)

    k2[0] = RHS_HE(rad+0.5*dr,p+0.5*k1[0],rho_plus_half_k1,m+0.5*k1[1])*dr
    k2[1] = RHS_MC(rad+0.5*dr,p+0.5*k1[0],rho_plus_half_k1,m+0.5*k1[1])*dr

    pnew = p + k2[0]
    mnew = m + k2[1]

    return (pnew,mnew)

def tov_integrate_RK3(rad,dr,p,rho,m):

    k1 = np.zeros(2)

    k1[0] = RHS_HE(rad,p,rho,m)*dr
    k1[1] = RHS_MC(rad,p,rho,m)*dr

    rho_plus_half_k1 = press_to_rho(p+0.5*k1[0])

    k2 = np.zeros(2)

    k2[0] = RHS_HE(rad+0.5*dr,p+0.5*k1[0],rho_plus_half_k1,m+0.5*k1[1])*dr
    k2[1] = RHS_MC(rad+0.5*dr,p+0.5*k1[0],rho_plus_half_k1,m+0.5*k1[1])*dr

    rho_minus_k1_plus_2k2 = press_to_rho(p-k1[0]+2*k2[0])

    k3 = np.zeros(3)

    k3[0] = RHS_HE(rad+dr,p-k1[0]+2*k2[0],rho_minus_k1_plus_2k2,m-k1[1]+2*k2[1])*dr
    k3[1] = RHS_MC(rad+dr,p-k1[0]+2*k2[0],rho_minus_k1_plus_2k2,m-k1[1]+2*k2[1])*dr

    pnew = p + (1.0/6.0)*(k1[0]+4*k2[0]+k3[0])
    mnew = m + (1.0/6.0)*(k1[1]+4*k2[1]+k3[1])

    return (pnew,mnew)

def tov_integrate_RK4(rad,dr,p,rho,m):

    k1 = np.zeros(2)

    k1[0] = RHS_HE(rad,p,rho,m)*dr
    k1[1] = RHS_MC(rad,p,rho,m)*dr

    rho_plus_half_k1 = press_to_rho(p+0.5*k1[0])

    k2 = np.zeros(2)

    k2[0] = RHS_HE(rad+0.5*dr,p+0.5*k1[0],rho_plus_half_k1,m+0.5*k1[1])*dr
    k2[1] = RHS_MC(rad+0.5*dr,p+0.5*k1[0],rho_plus_half_k1,m+0.5*k1[1])*dr

    rho_plus_half_k2 = press_to_rho(p+0.5*k2[0])

    k3 = np.zeros(2)

    k3[0] = RHS_HE(rad+0.5*dr,p+0.5*k2[0],rho_plus_half_k2,m+0.5*k2[1])*dr
    k3[1] = RHS_MC(rad+0.5*dr,p+0.5*k2[0],rho_plus_half_k2,m+0.5*k2[1])*dr

    rho_plus_k3 = press_to_rho(p+k3[0])

    k4 = np.zeros(2)

    k4[0] = RHS_HE(rad+dr,p+k3[0],rho_plus_k3,m+k3[1])*dr
    k4[1] = RHS_MC(rad+dr,p+k3[0],rho_plus_k3,m+k3[1])*dr

    pnew = p + (1.0/6.0)*(k1[0]+2*k2[0]+2*k3[0]+k4[0])
    mnew = m + (1.0/6.0)*(k1[1]+2*k2[1]+2*k3[1]+k4[1])

    return (pnew,mnew)

#######################################

# set up grid
npoints = 10000
radmax = 2.0e8 # 2000 km
radius = np.linspace(1.0e-8,radmax,npoints)
dr = radius[1]-radius[0]

# set up variables
press = np.zeros(npoints)
rho   = np.zeros(npoints)
mass  = np.zeros(npoints)

# set up central values
rho[0]   = 1.0e10
press[0] = polyK * rho[0]**polyG
mass[0]  = 0.0

# set up termination criterion
press_min = 1.0e-10 * press[0]

nsurf = 0
for n in range(npoints-1):
    
    #(press[n+1],mass[n+1]) = tov_integrate_FE(radius[n],dr,press[n],rho[n],mass[n])
    #(press[n+1],mass[n+1]) = tov_integrate_RK2(radius[n],dr,press[n],rho[n],mass[n])
    #(press[n+1],mass[n+1]) = tov_integrate_RK3(radius[n],dr,press[n],rho[n],mass[n])
    (press[n+1],mass[n+1]) = tov_integrate_RK4(radius[n],dr,press[n],rho[n],mass[n])

    # check for termination criterion
    if(press[n+1] < press_min and nsurf==0):
        nsurf = n

    if(n+1 > nsurf and nsurf > 0):
        press[n+1] = press[nsurf]
        rho[n+1]   = rho[nsurf]
        mass[n+1]  = mass[nsurf]

    # invert the EOS to get density
    rho[n+1] = press_to_rho(press[n+1])


print radius[nsurf]/1.0e5
print mass[nsurf]/msun


radius /= 1.0e5
mass /= msun
rho /= 1.0e10


import matplotlib.pyplot as plt
from matplotlib import rc

# set the axes and ticks of the plots
rc('axes',linewidth=2,labelsize=18);
rc('lines',linewidth=2);
rc('xtick',labelsize=18);
rc('ytick',labelsize=18);
rc('xtick.major',size=10,pad=8,width=1.5);
rc('xtick.minor',size=5);
rc('ytick.major',size=10,pad=8,width=1.5);
rc('ytick.minor',size=5);
rc('legend',fontsize=18);
rc('text', usetex=True)

# set plot position
myfig=plt.figure(0)
myfig.subplots_adjust(left=0.12)
myfig.subplots_adjust(bottom=0.14)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.88)

ax1 = myfig.add_subplot(111)
p1, = ax1.plot(radius, rho, 'b-')
ax1.set_xlabel('Radius (km)')
ax1.set_ylabel(r'Density ($\rm 10^{10}~g~cm^{-3}$)')

ax2 = ax1.twinx()
p2, = ax2.plot(radius, mass, 'g-')
ax2.set_ylabel('Mass ($M_{\odot}$)')
plt.legend([p1,p2,],["Density","Mass"],loc=4)
plt.show()
