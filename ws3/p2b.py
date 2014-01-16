import numpy as np
from scipy import special as sp

def g(x):
    return (x**2)/(1.0+np.exp(x/20.0))

energy=np.arange(0,155,5)
interval=len(energy)-1
integrate=np.zeros(interval)

for i in np.arange(interval):
    a=energy[i]; b=energy[i+1];
    roots,weights=sp.p_roots(10);
    x=a+0.5*(b-a)*(roots+1);
    integrate[i]=np.sum(weights*0.5*(b-a)*g(x))

print "i","E_i","n_i","n_i/dE"
for i in np.arange(interval):
    print i+1,energy[i],integrate[i],integrate[i]/5.
