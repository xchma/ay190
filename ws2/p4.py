import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
from matplotlib import rc

data=np.loadtxt("tb2.txt");
time=data[:,0]; mag=data[:,1];

# Lagrange interpolation
lagrange=interp.lagrange(time,mag)
t=np.linspace(0.0,1.0,200)
lag_interp=lagrange(t)

# piecewise linear
linear=interp.interp1d(time,mag,kind='linear')
linear_interp=linear(t)

# piecewise cubic
quad=interp.interp1d(time,mag,kind='quadratic')
quad_interp=quad(t)

# cubic Hermite interpolation
cubic=interp.interp1d(time,mag,kind='cubic')
cubic_interp=cubic(t)

# natural cubic spline
tck=interp.splrep(time,mag,k=3)
cspline_interp=interp.splev(t,tck)

# set the axes and ticks of the plots
rc('axes',linewidth=2,labelsize=16);
rc('xtick',labelsize=16);
rc('ytick',labelsize=16);
rc('xtick.major',size=10,pad=8,width=1.5);
rc('ytick.major',size=10,pad=8,width=1.5);
rc('legend',fontsize=16);

# set plot position
myfig=plt.figure(0)
myfig.subplots_adjust(left=0.12)
myfig.subplots_adjust(bottom=0.14)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.97)

plt.plot(time,mag,'ko',markersize=8.0)
plt.plot(t,lag_interp,'b-',t,linear_interp,'g-',t,quad_interp,'y', \
	t,cubic_interp,'r-',t,cspline_interp,'c-',linewidth=2)
plt.legend(["data","Lagrange","Linear","Quadratic","Cubic Hermite","Cubic Spline"],loc='best')
plt.xlabel("Time [days]");
plt.ylabel("Apparent Magnitude")
plt.savefig("fig3.pdf")


myfig=plt.figure(1)
myfig.subplots_adjust(left=0.12)
myfig.subplots_adjust(bottom=0.14)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.97)

plt.plot(time,mag,'ko',markersize=8.0)
plt.plot(t,lag_interp,'b-',t,linear_interp,'g-',t,quad_interp,'y', \
        t,cubic_interp,'r-',t,cspline_interp,'c-',linewidth=2)
plt.legend(["data","Lagrange","Linear","Quadratic","Cubic Hermite","Cubic Spline"],loc='best')
plt.xlabel("Time [days]")
plt.ylim([0.0,0.7])
plt.savefig("fig4.pdf")
