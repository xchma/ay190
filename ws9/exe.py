import numpy as np
import time
import scipy.linalg as linalg
import scipy.sparse.linalg as lin

import gaussian

def eqn(i,save_flag=1):

    # define file names
    A_file="LSE"+str(i)+"_m.dat";
    b_file="LSE"+str(i)+"_bvec.dat";
    x_file="LSE"+str(i)+"_soln.dat";

    # load files
    A=np.loadtxt(A_file);
    b=np.loadtxt(b_file);

    # gauss elimination
    start_time=time.time();
    soln=gaussian.gaussian(A,b);
    finish_time=time.time();

    # write to file
    if(save_flag==1): np.savetxt(x_file,soln);
    return finish_time-start_time


# for i in np.arange(1,6,1):
#   print i, eqn(i)

def eqn(i):
    
    # define file names
    A_file="LSE"+str(i)+"_m.dat";
    b_file="LSE"+str(i)+"_bvec.dat";
    x_file="LSE"+str(i)+"_soln.dat";
    
    # load files
    A=np.loadtxt(A_file);
    b=np.loadtxt(b_file);

    start_time=time.time();
    lu_and_piv=linalg.lu_factor(A);
    result=linalg.lu_solve(lu_and_piv,b);
    finish_time=time.time();    
    print finish_time-start_time

