import numpy as np
import matplotlib.pyplot as plt

x=np.zeros(16); x=np.float32(x);
x[0]=np.float32(1.0); x[1]=np.float32(1.0/3.0);

coeff1=np.float32(13.0)/np.float32(3.0);
coeff2=np.float32(4.0)/np.float32(3.0);

for n in np.arange(2,16,1):
    x[n]=coeff1*x[n-1]-coeff2*x[n-2];

xx=(1.0/3.0)**np.arange(0,16,1);
for n in np.arange(0,16,1):
    print n,x[n],xx[n],x[n]-xx[n],(x[n]-xx[n])/xx[n]
