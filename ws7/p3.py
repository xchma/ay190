import numpy as np
import numpy.random as random
import matplotlib.pyplot as plt
from matplotlib import rc

def f(x):
    return x**2+1.0

def MC_integrate(f,N,xmin,xmax):
    x=xmin+(xmax-xmin)*random.rand(N)
    y=f(x)
    return np.sum(y)*(xmax-xmin)/N

N_list=np.array([])
error_list=np.array([])
N=1000; xmin=2; xmax=3; exact=22.0/3.0
while (N<1.0e8):
    N_list=np.append(N_list,N)
    MC=MC_integrate(f,N,xmin,xmax)
    error=MC-exact
    error_list=np.append(error_list,np.absolute(error))
    print N,MC,error
    N*=2

# set the axes and ticks of the plots
rc('axes',linewidth=2,labelsize=20);
rc('lines',linewidth=2);
rc('xtick',labelsize=20);
rc('ytick',labelsize=20);
rc('xtick.major',size=10,pad=8,width=1.5);
rc('xtick.minor',size=5);
rc('ytick.major',size=10,pad=8,width=1.5);
rc('ytick.minor',size=5);
rc('legend',fontsize=16);
rc('text', usetex=True)

# set plot position
myfig=plt.figure(0)
myfig.subplots_adjust(left=0.12)
myfig.subplots_adjust(bottom=0.14)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.97)

plt.plot(N_list**(0.5),error_list,'bo-')
plt.xlabel("$\sqrt{N}$")
plt.ylabel("Error")
plt.show()