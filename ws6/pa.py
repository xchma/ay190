import numpy as np
from dft import *

def f(x):
    coe=1/np.sqrt(2*np.pi);
    return coe*np.exp(-0.5*x**2)

x=np.linspace(-2,2,30)
fx=f(x);
y1=dft(fx); y2=np.fft.fft(fx)

for i in np.arange(len(x)):
    if (y1[i].imag<0):
    	print '%.9f%.9f$i$' % (y1[i].real, y1[i].imag), \
	    '%.9f%.9f$i$' % (y2[i].real, y2[i].imag)
    if (y1[i].imag>=0):
        print '%.9f+%.9f$i$' % (y1[i].real, y1[i].imag), \
            '%.9f+%.9f$i$' % (y2[i].real, y2[i].imag)

