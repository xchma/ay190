import numpy as np
import numpy.random as random

def single_exp(N):
    birth=random.randint(365,size=N)
    birth_unique=np.unique(birth)
    return (len(birth_unique)<len(birth))

def multi_exp(N,N_exp):
    result=np.array([]);
    for i in np.arange(N_exp):
	result=np.append(result,single_exp(N));
    N_true=len(result[result==True]);
    return np.float(N_true)/np.float(N_exp)

N_list=np.arange(21,26,1);
prob_list=np.array([]);
for N in N_list:
    prob=multi_exp(N,20000)
    prob_list=np.append(prob_list,prob)
    print N,prob
