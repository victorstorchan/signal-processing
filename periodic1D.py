import numpy as np
import matplotlib.pyplot as plt
from cmath import polar
from math import sqrt
#definition of the boxcars signals

def boxcar(x,i):
    if x-i<-1 or x-i>1:
        return 0
    else:
        return 1
x= np.arange(-2.,2.,0.05)
n=len(x)
print(n)
True_signal=np.zeros(n)
for i in range(n):
    True_signal[i]=boxcar(x[i],0)
#plt.plot(x,True_signal)
#plt.axis([-2,2,-1,2])
#plt.show()
#definitions of the shifted signals
y=np.zeros(n,dtype=complex)
y2=np.zeros(n,dtype=complex)
base=np.zeros(n,dtype=complex)
vector_of_shift=[0,3,10,30]#shifts are integer in discrete version
len_shift=len(vector_of_shift)
#signal with shift:
shifted_signals=np.zeros((len_shift,n),dtype=complex)
shifted_signals_1=np.zeros((len_shift,n),dtype=complex)
for k in range(n):
    base[k]=boxcar(x[k],0)
max_shift=max(vector_of_shift)
base_period=np.lib.pad(base, (max_shift, 0), 'wrap')
for s in range(len_shift):
    for k in range(n):
        if k-vector_of_shift[s]<0:
            y[k]=base_period[max_shift-vector_of_shift[s]-1+k]
            y2[k]=base_period[max_shift-vector_of_shift[s]-1+k]*np.exp(2J*np.pi*k/n)
        else:
            y[k]=boxcar(x[k-vector_of_shift[s]],0)
            y2[k]=boxcar(x[k-vector_of_shift[s]],0)*np.exp(2J*np.pi*k/n)
    
    randvect=np.random.normal(0,0.1,n)
    shifted_signals[s] =y#+ randvect
    shifted_signals_1[s]=y2#+ randvect


A=np.fft.fft(shifted_signals)
A_1=np.fft.fft(shifted_signals_1).conjugate()
A_star=np.zeros((len_shift,n),dtype=complex)
for i in range(len_shift):
    A_star[i] = A[i]*A_1[i]
    
A_star_matrix=np.matrix(A_star)
A_star_transpose=A_star_matrix.getH()
A_prod1=A_star_matrix*A_star_transpose
A_prod=A_prod1/A_prod1[0,0]

(V,sigma,V_star)=np.linalg.svd(A_prod,full_matrices=1)
v1=V_star[0].getH()

#the shifts are recovered:
output=np.zeros(len_shift,dtype=complex)
for i in range(len_shift):
    output[i]=-n*polar(-v1[i,0])[1]/(2*np.pi)
output
