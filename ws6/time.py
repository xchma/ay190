import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from timeit import timeit

setup = \
 "import numpy as np; import dft; x=np.random.random(%d)" 
setup_fft = \
 "import numpy as np; x=np.random.random(%d)"

N_list=np.arange(10,110,10);
dft_time=np.zeros(len(N_list));
fft_time=np.zeros(len(N_list));
for i in np.arange(len(N_list)):
    dft_time[i]=timeit("dft.dft(x)",number=100,setup=setup % N_list[i])
    fft_time[i]=timeit("np.fft.fft(x)",number=100, 
	setup=setup_fft % N_list[i])

# set the axes and ticks of the plots
rc('axes',linewidth=2,labelsize=16);
rc('xtick',labelsize=14);
rc('ytick',labelsize=14);
rc('xtick.major',size=10,pad=8,width=1.5);
rc('ytick.major',size=10,pad=8,width=1.5);
rc('legend',fontsize=16);

myfig=plt.figure(0)
myfig.subplots_adjust(left=0.13)
myfig.subplots_adjust(bottom=0.14)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

plt.plot(N_list**2,dft_time,"bo-",markersize=8,linewidth=2)
plt.xlabel(r"$N^2$")
plt.ylabel("time")
plt.savefig("fig1.pdf",format='pdf')

myfig=plt.figure(1)
myfig.subplots_adjust(left=0.13)
myfig.subplots_adjust(bottom=0.14)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

plt.plot(N_list*(np.log2(N_list)),fft_time,"bo-",markersize=8,linewidth=2)
plt.xlabel(r"$N\log_2(N)$")
plt.savefig("fig2.pdf",format='pdf')