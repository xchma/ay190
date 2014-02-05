import numpy as np

def gaussian(A,b):

    # test the input arrays
    A_shape=A.shape; b_shape=b.shape;
    if ((len(A_shape)!=2) | (len(b_shape)!=1) | (A_shape[0]!=A_shape[1]) | (A_shape[0]!=b_shape[0])):
        print "Input Error";
	return

    # begin the loop of elimination
    N=A_shape[0]
    for i in np.arange(N):
        # find the max line
        which_line=i; pivot=np.absolute(A[i,i]);
        for j in np.arange(i+1,N,1):
            if (np.absolute(A[j,i])>pivot):
                pivot=np.absolute(A[j,i]);
                which_line=j;
        
        # switch line
        temp_line=np.copy(A[which_line]); temp=b[which_line];
        A[which_line]=A[i]; b[which_line]=b[i];
        A[i]=temp_line; b[i]=temp;

        # elimination
        conv=A[i,i];
        if (conv==0):
            print "This system is not solvable"
            return
        A[i]=A[i]/conv; b[i]=b[i]/conv;

        for j in np.arange(i+1,N,1):
            conv=A[j,i]/A[i,i];
            A[j]=A[j]-conv*A[i];
            b[j]=b[j]-conv*b[i];

    # begin the trace back
    for i in np.arange(N):
        for j in np.arange(i+1,N,1):
            conv=A[N-1-j,N-1-i];
            A[N-1-j]=A[N-1-j]-conv*A[N-1-i];
            b[N-1-j]=b[N-1-j]-conv*b[N-1-i];

    return b

