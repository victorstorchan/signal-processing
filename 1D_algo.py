# coding: utf-8


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

True_signal=np.zeros(n)
for i in range(n):
    True_signal[i]=boxcar(x[i],0)
#plt.plot(x,True_signal)
#plt.axis([-2,2,-1,2])
#plt.show()
#definitions of the shifted signals
y=np.zeros(n,dtype=complex)
y2=np.zeros(n,dtype=complex)

vector_of_shift=[0,2,4,10]#shifts are integer in discrete version
len_shift=len(vector_of_shift)
#signal with shift:
shifted_signals=np.zeros((len_shift,n),dtype=complex)
shifted_signals_1=np.zeros((len_shift,n),dtype=complex)

for s in range(len_shift):
    for k in range(n):
        if k<vector_of_shift[s]:
            y[k]=0
            y2[k]=0
        else:
            y[k]=boxcar(x[k-vector_of_shift[s]],0)
            y2[k]=boxcar(x[k-vector_of_shift[s]],0)*np.exp(2J*np.pi*k/n)
    randvect=np.random.normal(0,0.1,n)
    shifted_signals[s] =y+ randvect
    shifted_signals_1[s]=y2+ randvect
y_noisy_1=shifted_signals[0]
y_noisy_2=shifted_signals[1]
y_noisy_3=shifted_signals[2]
y_noisy_4=shifted_signals[3]
#plot of the signals:
#plt.plot(x,y_noisy_3)
#plt.axis([-2,2,-1,2])
#plt.show()
# Four axes, returned as a 2-d array
f, axarr = plt.subplots(2, 2)
axarr[0, 0].plot(x, y_noisy_1)
axarr[0, 0].set_title('shift_val = 0')
axarr[0, 1].plot(x, y_noisy_2)
axarr[0, 1].set_title('shift_val = 2')
axarr[1, 0].plot(x, y_noisy_3)
axarr[1, 0].set_title('shift_val = 4')
axarr[1, 1].plot(x, y_noisy_4)
axarr[1, 1].set_title('shift_val = 10')
plt.show()



A=np.fft.fft(shifted_signals)
A_1=np.fft.fft(shifted_signals_1).conjugate()
A_star=np.zeros((len_shift,n),dtype=complex)
for i in range(len_shift):
    A_star[i] = A[i]*A_1[i]
    
A_star_matrix=np.matrix(A_star)
A_star_transpose=A_star_matrix.getH()
A_prod1=A_star_matrix*A_star_transpose
A_prod=A_prod1/A_prod1[0,0]

(V,sigma,V_star)=np.linalg.svd(A_prod,full_matrices=0)
v1=V_star[0].getH()




#the shifts are recovered:
output=np.zeros(len_shift,dtype=complex)
for i in range(len_shift):
    output[i]=-n*polar(-v1[i,0])[1]/(2*np.pi)
output



#The signal is recovered:
recover_A=np.zeros((len_shift,n),dtype=complex)
for i in range(len_shift):
    for j in range(n):
        recover_A[i,j]=A[i,j]/np.exp(-2J*np.pi*j*output[i]/n)



A_final=np.fft.ifft(recover_A)



#the more signals we have, the less the noise appears at the end

y_recov=np.zeros(n,dtype=complex)
for i in range(len_shift):
    y_recov+=A_final[i]
y_recov=y_recov/(len_shift)



#plot of the signals:
plt.plot(x,y_recov)
plt.plot(x,True_signal,'r')
plt.axis([-2,2,-1,2])
plt.show()

