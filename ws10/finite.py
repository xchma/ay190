import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse.linalg as splin

def finite_M(x,g,p,q):
    
    npoints=len(x)-2
    h=x[1]-x[0] # only valid for uniform grid
    
    A=np.zeros((npoints,npoints))
    for i in np.arange(npoints):
        A[i,i]=-2+h**2*q(x[i+1])
    for i in np.arange(npoints-1):
        A[i,i+1]=1+0.5*h*p(x[i+1])
    for i in np.arange(1,npoints,1):
        A[i,i-1]=1-0.5*h*p(x[i+1])

    return A

def finite_b(x,A,B,g,p):

    npoints=len(x)-2
    h=x[1]-x[0] # only valid for uniform grid

    b=np.zeros(npoints)
    for i in np.arange(npoints):
        b[i]=h**2*g(x[i+1])

    b[0]=b[0]-A*(1-0.5*h*p(x[1]))
    b[-1]=b[-1]-B*(1+0.5*h*p(x[npoints]))
                   
    return b

def g(x):
    return 12*x-4

def p(x):
    return 0*x

def q(x):
    return 0*x

xmin=0.0
xmax=1.0
A=0.0
B=0.1

def test(npoints):
    x=np.linspace(xmin,xmax,npoints+1)
    y=np.zeros(len(x))
    y[0]=A
    y[-1]=B

    M=finite_M(x,g,p,q)
    b=finite_b(x,A,B,g,p)

    y[1:-1]=splin.spsolve(M,b)

    return x,y

x1,y1=test(6)
x2,y2=test(10)

xinf=np.linspace(xmin,xmax,1000)

plt.plot(x1,y1,'b-',x2,y2,'g-')
plt.plot(xinf,2.0*xinf**3-2*xinf**2+0.1*xinf,"r-")
plt.xlabel('x')
plt.ylabel('y')
plt.title('Finite Difference Method')
plt.legend(["6 points","10 points","true solution"],loc='best')
plt.show()

