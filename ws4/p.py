import numpy as np

# phi=omega*t, n are parameters
def f(x,phi,e):
    return x-phi-e*np.sin(x)

# need two initial x-values
def secant(f,x1,x2):
    delta=x2-x1; count=0;
    while (np.absolute(delta)>=1.0e-10):
	x3=x2-f(x2)*delta/(f(x2)-f(x1));
	x1=x2; x2=x3; delta=x2-x1;
	count+=1;
    return x2,count

# star to calculate
# define some parameters
for e in [0.0167,0.99999]:
    a=1.0; b=np.sqrt(a**2*(1-e**2)); # semi-axis
    T=365.25635; # period
    omega=2*np.pi/T

    for t in np.array([91.0,182.0,273.0]):
	phi=omega*t;

	def g(x): 
	    return f(x,phi,e)

	E,count=secant(g,phi,phi+1.0e-3)
	x=a*np.cos(E); y=b*np.sin(E);

	print t,E,x,y,count
