
# coding: utf-8



#2-D signals
import numpy as np
import matplotlib.pyplot as plt
from cmath import polar
from PIL import Image 
from random import randint




def create_shift(K):
    l=[]
    for i in range(K):
        l.append((randint(0,40),randint(0,40)))
    return l




#K=number of shift, sigma=noise magnitude
def signal_to_noise(K,noise_level):
    image_file = Image.open("/Users/victorstorchan/Desktop/RA/ra/apple_original signal.png") # open colour image
    image_signal = image_file.convert('L') # convert image to black and white
    signal=np.asarray(image_signal)
    (a,b)=signal.shape
    randommat=np.zeros((400,400),dtype=complex)
    for i in range(400):
        for j in range(400):
            randommat[i][j]+=complex(np.random.normal(0,noise_level,1))
    signal_noisy=np.zeros((400,400),dtype=complex)
    for i in range(a):
        for j in range(b):
            signal_noisy[10+i][10+j]=signal[i][j]
    for i in range(400):
        for j in range(400):
            signal_noisy[i][j]+=randommat[i][j]
    m=signal_noisy.shape[0]
    n=signal_noisy.shape[1]
    vect_of_shift=create_shift(K)
    len_shift=len(vect_of_shift)
     #signal with shift:
    shifted_signals=[]
    shifted_signals_1=[]
    shifted_signals_2=[]
    for s in range(len_shift):
        y=np.zeros((n,m),dtype=complex)
        y1=np.zeros((n,m),dtype=complex)
        y2=np.zeros((n,m),dtype=complex)
        for k in range(m):
            for l in range(n):
                if (l<10+vect_of_shift[s][0] or l>315+vect_of_shift[s][0]) and (k<10+vect_of_shift[s][1] or k>324+vect_of_shift[s][1]):
                    y[k][l]=randommat[k][l]
                    y1[k][l]=randommat[k][l]*np.exp(2J*np.pi*k/m)
                    y2[k][l]=randommat[k][l]*np.exp(2J*np.pi*l/n)
                else:
                    y[k][l]=signal_noisy[k-vect_of_shift[s][0]][l-vect_of_shift[s][1]]
                    y1[k][l]=signal_noisy[k-vect_of_shift[s][0]][l-vect_of_shift[s][1]]*np.exp(2J*np.pi*k/m)
                    y2[k][l]=signal_noisy[k-vect_of_shift[s][0]][l-vect_of_shift[s][1]]*np.exp(2J*np.pi*l/n)
        shifted_signals.append(y)
        shifted_signals_1.append(y1)
        shifted_signals_2.append(y2)
    A=[]
    A_1=[]
    A_2=[]
    for i in range(len(shifted_signals)):
        A.append(np.matrix(np.fft.fft2(shifted_signals[i])))
        A_1.append(np.matrix(np.fft.fft2(shifted_signals_1[i])).conjugate())
        A_2.append(np.matrix(np.fft.fft2(shifted_signals_2[i])).conjugate())
    G1=[]
    G2=[]
    for s in range(len_shift):
        G1.append(np.multiply(A[s],A_1[s]))
        G2.append(np.multiply(A[s],A_2[s]))
    A_mat1=np.zeros((len_shift,n**2),dtype=complex)
    A_mat2=np.zeros((len_shift,n**2),dtype=complex)
    for s in range(len_shift):
        for i in range(n):
            for k in range(n):
                A_mat1[s][400*i+k]=G1[s][i,k]
                A_mat2[s][400*i+k]=G2[s][i,k]
    A_mat1_mat=np.matrix(A_mat1)
    A_mat2_mat=np.matrix(A_mat2)
    A_mat1_transpose=A_mat1_mat.getH()
    A_mat2_transpose=A_mat2_mat.getH()
    A_prod1=A_mat1_mat*A_mat1_transpose
    A_prod2=A_mat2_mat*A_mat2_transpose
    A_final1=A_prod1/A_prod1[0,0]
    A_final2=A_prod2/A_prod2[0,0]
    (V1,sigma1,V_star1)=np.linalg.svd(A_final1)
    (V2,sigma2,V_star2)=np.linalg.svd(A_final2)
    v1=V_star1[0].getH()
    v2=V_star2[0].getH()
    #the shifts are recovered:
    output1=np.zeros(len_shift)
    output2=np.zeros(len_shift)
    for i in range(len_shift):
        output1[i]=-n*polar(-v1[i,0])[1]/(2*np.pi).real
        output2[i]=-n*polar(-v2[i,0])[1]/(2*np.pi).real
    recover_A=[]
    for i in range(len_shift):
        M=np.matrix(np.zeros((m,n),dtype=complex))
        for l in range(m):
            for k in range(n):
                M[l,k]=A[i][l,k]/np.exp(-2J*np.pi*(l*output1[i]+k*output2[i])/n)
        recover_A.append(M)
    A_final=[]
    for i in range(len_shift):
        A_final.append(np.fft.ifft2(recover_A[i]))
    recov_signal=np.zeros((m,n))
    for i in range(m):
        for j in range(n):
            k=0
            for s in range(len_shift):
                k+=A_final[s][i][j].real
            recov_signal[i][j]=k/len_shift
    recov_signal1=recov_signal.astype(np.uint8)
    return np.linalg.norm(recov_signal1-im3)



#k is just the count of the number of iteration to see how the program is running, because 
#it is slow. 
result_matrix=np.zeros((100,100))
k=0
for i in range(2,100):
    for j in range(1,100):
        result_matrix[99-i][j]=signal_to_noise(i,j)
        k+=1
        print(k)
