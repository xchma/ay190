import numpy as np
import numpy.random as random
import matplotlib.pyplot as plt
from matplotlib import rc

def calc_pi(N):
    x=random.rand(N); y=random.rand(N);
    r=np.sqrt((x-0.5)**2+(y-0.5)**2);
    num_in_circle=len(r[r<0.5]);
    return 4*np.float(num_in_circle)/np.float(N)

err=np.array([]); 
N_list=np.array([]);
N=1000;
while (N<5.0e7):
    N_list=np.append(N_list,N);
    pi=calc_pi(N);
    err=np.append(err,np.absolute(pi-np.pi));
    print N,pi,pi-np.pi
    N *= 2;

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

plt.plot(N_list**(0.5),err,'bo-')
plt.xlabel("$\sqrt{N}$")
plt.ylabel("Error")
plt.show()
