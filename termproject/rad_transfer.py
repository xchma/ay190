import numpy as np
import ctypes

# define the grid
N_r = 20
N_theta = 20
r_min = 2.0
r_max = 4
theta_min = 0.0
theta_max = 1.6
N_ray = 20000

# assign some density
density = np.zeros((N_r, N_theta),dtype='f')
for i in np.arange(N_r):
    for j in np.arange(N_theta):
        density[i,j] = 1.0

# convert variables to ctypes
N_r_ctype = ctypes.c_int(N_r)
N_theta_ctype = ctypes.c_int(N_theta)
r_min_ctype = ctypes.c_float(r_min)
r_max_ctype = ctypes.c_float(r_max)
theta_min_ctype = ctypes.c_float(theta_min)
theta_max_ctype = ctypes.c_float(theta_max)
N_ray_ctype = ctypes.c_int(N_ray)
density_ctype = density.ctypes.data_as(ctypes.POINTER(ctypes.c_float))

# initialize the output pointer
r_final = (ctypes.c_float*N_ray)()
theta_final = (ctypes.c_float*N_ray)()

# find the C routine
exec_call = "./ray_trace/ray_trace.so"
ray_routine = ctypes.cdll[exec_call]

# call the routine
ray_routine.ray(N_r_ctype, N_theta_ctype, r_min_ctype, r_max_ctype, \
                theta_min_ctype, theta_max_ctype, \
                N_ray_ctype, density_ctype, \
                r_final, theta_final)

# get back the final coordinates
r_out = np.ctypeslib.as_array(r_final)
theta_out = np.ctypeslib.as_array(theta_final)

print r_out, theta_out