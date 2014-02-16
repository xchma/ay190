import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

# constants
G=6.67e-8

def euler(r,rho):

    N_grid=len(r)

    phi=np.zeros(N_grid)
    z=np.zeros(N_grid)
    M=np.zeros(N_grid)

    # inner boundary conditions, to avoid divergence
    phi[0]=1.0e-10 # this doesn't matter
    # in the case the grid starts from origin
    z[0]=1.0e-10
    M[0]=1.0e-10
    # if the grid starts from some finite radius
    if (r[0]>0):
        M[0]=(4.0/3.0)*np.pi*(r[0]**3)*rho[0]
        z[0]=G*M[0]/(r[0]**2)

    for i in np.arange(1,N_grid,1):
        dr=r[i]-r[i-1] # adaptive grid
        phi[i]=phi[i-1]+z[i-1]*dr
        z[i]=z[i-1]+(4*np.pi*G*rho[i-1]-2*z[i-1]/r[i-1])*dr
        M[i]=M[i-1]+4*np.pi*rho[i-1]*(r[i-1]**2)*dr

    # potential calibration
    phi_outer=-G*M[N_grid-1]/r[N_grid-1] # at the outer boundary
    phi=phi+phi_outer-phi[N_grid-1] # shifting

    return phi,z,M



# set the axes and ticks of the plots
rc('axes',linewidth=2,labelsize=20);
rc('lines',linewidth=2);
rc('xtick',labelsize=20);
rc('ytick',labelsize=20);
rc('xtick.major',size=10,pad=8,width=1.5);
rc('xtick.minor',size=5);
rc('ytick.major',size=10,pad=8,width=1.5);
rc('ytick.minor',size=5);
rc('legend',fontsize=20);
rc('text', usetex=True)

def plot_density_profile(r,rho):
    myfig=plt.figure(0)
    myfig.subplots_adjust(left=0.15)
    myfig.subplots_adjust(bottom=0.15)
    myfig.subplots_adjust(top=0.95)
    myfig.subplots_adjust(right=0.95)
    
    plt.plot(r,rho,'-')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([6.0e6,1.0e9])
    plt.ylim([1.0e4,2.0e10])
    plt.xlabel('radius (cm)')
    plt.ylabel(r'density ($\rm g~cm^{-3}$)')
    plt.savefig("fig1.pdf")

def test_euler():
    rout=1.0e5
    r=np.linspace(1.0e2,rout,200)
    rho=np.ones(len(r))
    
    phi,z,M=euler(r,rho)
    phi_ana=(2.0/3.0)*np.pi*G*(r**2-3*rout**2)
    
    myfig=plt.figure(0)
    myfig.subplots_adjust(left=0.15)
    myfig.subplots_adjust(bottom=0.15)
    myfig.subplots_adjust(top=0.95)
    myfig.subplots_adjust(right=0.95)
    
    plt.plot(r,phi,'b-',r,phi_ana,'g-')
    plt.xlabel('radius (cm)')
    plt.ylabel(r'$\Phi(r)$ (arbitrary units)')
    plt.legend(["Numerical","Analytical"],loc='best')
    plt.savefig("fig2.pdf")

data=np.loadtxt("presupernova.dat",usecols=(2,4))
radius=data[:,0]
density=data[:,1]

r=np.arange(6.0e6,1.0e9,1.0e6)
rho=np.interp(r,radius,density)
phi,z,M=euler(r,rho)

test_euler()

myfig=plt.figure(1)
myfig.subplots_adjust(left=0.15)
myfig.subplots_adjust(bottom=0.15)
myfig.subplots_adjust(top=0.95)
myfig.subplots_adjust(right=0.95)

plt.plot(r,phi,'b-')
plt.xlabel('radius (cm)')
plt.ylabel(r'$\Phi(r)$ (arbitrary units)')
plt.savefig("fig3.pdf")
