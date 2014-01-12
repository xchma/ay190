import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

def f(x):
  return x**3-5.0*x**2+x

def df(x):
  return 3.0*x**2-10.0*x+1.0

h1=0.01; h2=0.005;
x1=np.arange(-2,6,h1,dtype='float64');
x2=np.arange(-2,6,h2,dtype='float64');
y1=f(x1); y2=f(x2);
dy1=df(x1); dy2=df(x2); # real answer

# initialization
l1=len(x1); l2=len(x2);
dy1_forward=np.zeros(l1); dy1_central=np.zeros(l1);
dy2_forward=np.zeros(l2); dy2_central=np.zeros(l2);

# forward differencing with different resolution
for i in np.arange(0,l1-1,1):
    dy1_forward[i]=(y1[i+1]-y1[i])/h1
for i in np.arange(0,l2-1,1):
    dy2_forward[i]=(y2[i+1]-y2[i])/h2
# central differencing with different resolution
for i in np.arange(1,l1-1,1):
    dy1_central[i]=(y1[i+1]-y1[i-1])/(2*h1)
for i in np.arange(1,l2-1,1):
    dy2_central[i]=(y2[i+1]-y2[i-1])/(2*h2)

# absolute errors
err1_forward=dy1_forward-dy1;
err2_forward=dy2_forward-dy2;
err1_central=dy1_central-dy1;
err2_central=dy2_central-dy2;

# set the axes and ticks of the plots
rc('axes',linewidth=2,labelsize=16);
rc('xtick',labelsize=16);
rc('ytick',labelsize=16);
rc('xtick.major',size=10,pad=8,width=1.5);
rc('ytick.major',size=10,pad=8,width=1.5);
rc('legend',fontsize=16);

myfig=plt.figure(0)
myfig.subplots_adjust(left=0.13)
myfig.subplots_adjust(bottom=0.14)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

p1,=plt.plot(x1[0:l1-1],err1_forward[0:l1-1],'b-',linewidth=2.5);
p2,=plt.plot(x2[0:l2-1],err2_forward[0:l2-1],'g-',linewidth=2.5);
plt.legend((p1,p2),("forward h1=0.01","forward h2=0.005"),loc=2);
plt.xlabel('x');
plt.ylabel('absolute error');
plt.savefig("fig1.pdf",format='pdf');


myfig=plt.figure(1)
myfig.subplots_adjust(left=0.13)
myfig.subplots_adjust(bottom=0.14)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

p3,=plt.plot(x1[1:l1-1],err1_central[1:l1-1],'b-',linewidth=2.5);
p4,=plt.plot(x2[1:l2-1],err2_central[1:l2-1],'g-',linewidth=2.5);
plt.xlabel('x');
plt.legend((p3,p4),("central h1=0.01","central h2=0.005"),loc=10)
plt.savefig("fig2.pdf",format='pdf')
