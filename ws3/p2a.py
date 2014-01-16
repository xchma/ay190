import numpy as np
from scipy import special as sp

def f(x):
    temp=np.exp(x)
    return (x**2)*temp/(temp+1.0)

def laguerre(n,f):
    roots,weights=sp.l_roots(n);
    return np.sum(weights*f(roots))

print "n","I"
for n in np.arange(6,20,1):
    print n,laguerre(n,f)
