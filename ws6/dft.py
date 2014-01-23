import numpy as np

def dft(x):
    N=len(x);
    A=np.zeros((N,N),dtype=complex);
    w=np.exp(-2*np.pi*(1j)/np.float(N))
    for k in np.arange(N):
	for j in np.arange(N):
	    A[k,j]=w**(k*j);
    return np.dot(A,x)
