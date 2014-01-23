import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from matplotlib import rc

def linear_fitting(x,y,sigma):
    S=np.sum(1/sigma**2);
    sx=np.sum(x/sigma**2);
    sy=np.sum(y/sigma**2);
    sxx=np.sum(x**2/sigma**2);
    sxy=np.sum(x*y/sigma**2);
    a1=(sy*sxx-sx*sxy)/(S*sxx-sx**2);
    a2=(S*sxy-sx*sy)/(S*sxx-sx**2);
    return a1,a2

# read the data
data=ascii.read("m_sigma_table.dat",readme="m_sigma_ReadMe.dat");
sigma=np.array(data["sigma*"]);
sigma_err=np.array(data["e_sigma*"]);
# convert velocity dispersion to logarithmic scale
sigma_err=sigma_err/(sigma*np.log(10));
sigma=np.log10(sigma);
mbh=np.array(data["logM"]);
mbh_error=np.array(data["e_logM"]);

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

# fitting without error
x=sigma; y=mbh; error=np.ones(len(x));
a1,a2=linear_fitting(x,y,error);
sigma_fitting=np.linspace(1.4,2.6,200);
mbh_fitting=a1+a2*sigma_fitting;

plt.plot(sigma,mbh,'bo',markersize=8.0);
plt.plot(sigma_fitting,mbh_fitting,'k--');

# fitting with error
x=sigma; y=mbh; error=mbh_error;
a1,a2=linear_fitting(x,y,error);
sigma_fitting=np.linspace(1.4,2.6,200);
mbh_fitting=a1+a2*sigma_fitting;

plt.plot(sigma_fitting,mbh_fitting,'k-');

mbh_fitting=-0.64+3.69*sigma_fitting;
plt.plot(sigma_fitting,mbh_fitting,'k-.');
plt.xlabel(r"$\log(\sigma_*/\rm km~s^{-1})$");
plt.ylabel(r"$\log(M_{\rm BH}/M_{\odot})$");
plt.legend(["data","fitting without error","fitting with error", \
	"Greene \& Ho (2006)"],loc='best');
plt.xlim([1.4,2.6]);
plt.savefig("fig.pdf",format='pdf')

