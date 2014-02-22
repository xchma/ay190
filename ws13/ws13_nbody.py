#!/usr/bin/env python

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as mpl3d

# global constants
ggrav = 6.67e-8
msun  = 1.99e33
seconds_per_year = 24.*3600*365 # roughly
cm_per_pc = 3.1e18
distance_to_sgrAstar = 8e3 * cm_per_pc

# system parameters
# run sun-earth system, using these parameters

initial_data_file = "sun_earth.asc.txt"
distance_unit_to_cm = 1.
time_unit_to_s = 1.
mass_unit_to_g = 1.
Nsteps=5000
t0 = 0
t1 = 5 * seconds_per_year

# run SgrA* system, using these parameters

#initial_data_file = "sgrAstar.asc.txt"
#distance_unit_to_cm = distance_to_sgrAstar * (1.0/3600.0) * (np.pi/180.0)
#time_unit_to_s = seconds_per_year
#mass_unit_to_g = msun
#Nsteps = 10000
#t0 = 0
#t1 = 100 * seconds_per_year

dt = (t1-t0)/Nsteps

final_data_file = "final_positions_sgrAstar.asc.txt"

def NbodyRHS(u,mass):
    (m,n) = u.shape
    rhs = np.zeros((m,n))

    rhs[:,0] = u[:,3]
    rhs[:,1] = u[:,4]
    rhs[:,2] = u[:,5]

    for i in np.arange(m):
        force = np.array([0.0,0.0,0.0])
        for j in np.arange(m):
            if (j==i): continue
            rel_pos = np.array([u[j,0]-u[i,0],u[j,1]-u[i,1],u[j,2]-u[i,2]])
            mol = np.sqrt(rel_pos[0]**2+rel_pos[1]**2+rel_pos[2]**2)
            force += ggrav*mass[j]*rel_pos/mol**3
        rhs[i,3:6]=force

    return rhs

def NbodyRK4(u,mass,time,dt):
    
    k1 = NbodyRHS(u,mass)*dt
    k2 = NbodyRHS(u+0.5*k1,mass)*dt
    k3 = NbodyRHS(u+0.5*k2,mass)*dt
    k4 = NbodyRHS(u+k3,mass)*dt

    u += (1.0/6.0)*(k1+2*k2+2*k3+k4)
    
    return u

def TotalEnergy(u,mass,time):
    (m,n) = u.shape

    energy = 0.0

    for i in np.arange(m):
        energy += 0.5*mass[i]*np.sqrt(u[i,3]**2+u[i,4]**2+u[i,5]**2)

    for i in np.arange(m):
        for j in np.arange(m):
            if (j==i): continue
            rel_pos = np.array([u[j,0]-u[i,0],u[j,1]-u[i,1],u[j,2]-u[i,2]])
            mol = np.sqrt(rel_pos[0]**2+rel_pos[1]**2+rel_pos[2]**2)
            energy += -ggrav*mass[i]*mass[j]/mol

    return energy


# main program
plt.ion()

(x,y,z,vx,vy,vz,mass) = np.loadtxt(initial_data_file, unpack = True)


# convert from unitis in initial data file to cgs
x *= distance_unit_to_cm
y *= distance_unit_to_cm
z *= distance_unit_to_cm
vx *= distance_unit_to_cm / time_unit_to_s
vy *= distance_unit_to_cm / time_unit_to_s
vz *= distance_unit_to_cm / time_unit_to_s
mass *= mass_unit_to_g

xmin = np.amin(x)
xmax = np.amax(x)
ymin = np.amin(y)
ymax = np.amax(y)
zmin = np.amin(z)
zmax = np.amax(z)
rmax = 2.5*max(abs(xmin),abs(xmax),abs(ymin),abs(ymax),abs(zmin),abs(zmax))

# use a single state vector to simplify the ODE code
# indices:
# u[:,0] = x
# u[:,1] = y
# u[:,2] = z
# u[:,3] = vx
# u[:,4] = vy
# u[:,5] = vz
u = np.array((x,y,z,vx,vy,vz)).transpose()

for it in range(0, Nsteps):
    time = t0 + it * dt
    u = NbodyRK4(u,mass,time,dt)
    if it % max(1,Nsteps/100) == 0:
      print "it = %d, time = %g years, energy = %g" % \
            (it, time / seconds_per_year,
             TotalEnergy(u,mass,time))
      plt.clf()
      fig = plt.gcf()
      ax = mpl3d.Axes3D(fig)
      ax.scatter(u[:,0],u[:,1],u[:,2])
      ax.set_xlim((-rmax,rmax))
      ax.set_ylim((-rmax,rmax))
      ax.set_zlim((-rmax,rmax))
      plt.draw()

# output result
file_header = "1:x 2:y 3:z 4:vx 5:vy 6:vz 7:mass"
np.savetxt(final_data_file, u, header=file_header)
