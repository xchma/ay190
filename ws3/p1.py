import numpy as np

def f1(x):
    return np.sin(x)

def f2(x):
    return x*np.sin(x)

def mid_point(x,f):
    intervals=len(x)-1;
    integrate=np.zeros(intervals);
    for i in np.arange(intervals):
	integrate[i]=(x[i+1]-x[i])*f((x[i+1]+x[i])*0.5)
    return np.sum(integrate)

def trapezoidal(x,f):
    intervals=len(x)-1;
    integrate=np.zeros(intervals);
    for i in np.arange(intervals):
        integrate[i]=0.5*(x[i+1]-x[i])*(f(x[i+1])+f(x[i]))
    return np.sum(integrate)

def simpson(x,f):
    intervals=len(x)-1;
    integrate=np.zeros(intervals);
    coeff=1.0/6.0;
    for i in np.arange(intervals):
        integrate[i]=coeff*(x[i+1]-x[i])*(f(x[i+1])+f(x[i])+4*f((x[i+1]+x[i])*0.5))
    return np.sum(integrate)

# generate the grid points
x1=np.linspace(0,np.pi,100);
x2=np.linspace(0,np.pi,200);
x3=np.linspace(0,np.pi,500);
x4=np.linspace(0,np.pi,1000);

# do the integration for part a
mp1=mid_point(x1,f1); mp2=mid_point(x2,f1);
mp3=mid_point(x3,f1); mp4=mid_point(x4,f1);
tr1=trapezoidal(x1,f1); tr2=trapezoidal(x2,f1);
tr3=trapezoidal(x3,f1); tr4=trapezoidal(x4,f1);
si1=simpson(x1,f1); si2=simpson(x2,f1);
si3=simpson(x3,f1); si4=simpson(x4,f1);

# print the results
print "Mid-point",mp1,mp2,mp3,mp4
print "Error",mp1-2.0,mp2-2.0,mp3-2.0,mp4-2.0
print "Trapezoidal",tr1,tr2,tr3,tr4
print "Error",tr1-2.0,tr2-2.0,tr3-2.0,tr4-2.0
print "Simpson",si1,si2,si3,si4
print "Error",si1-2.0,si2-2.0,si3-2.0,si4-2.0

# do the integration for part b
mp1=mid_point(x1,f2); mp2=mid_point(x2,f2);
mp3=mid_point(x3,f2); mp4=mid_point(x4,f2);
tr1=trapezoidal(x1,f2); tr2=trapezoidal(x2,f2);
tr3=trapezoidal(x3,f2); tr4=trapezoidal(x4,f2);
si1=simpson(x1,f2); si2=simpson(x2,f2);
si3=simpson(x3,f2); si4=simpson(x4,f2);

# print the results
print "Mid-point",mp1,mp2,mp3,mp4
print "Error",mp1-np.pi,mp2-np.pi,mp3-np.pi,mp4-np.pi
print "Trapezoidal",tr1,tr2,tr3,tr4
print "Error",tr1-np.pi,tr2-np.pi,tr3-np.pi,tr4-np.pi
print "Simpson",si1,si2,si3,si4
print "Error",si1-np.pi,si2-np.pi,si3-np.pi,si4-np.pi
