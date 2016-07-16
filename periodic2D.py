#2-D signals
import numpy as np
import matplotlib.pyplot as plt
from cmath import polar
from PIL import Image

def shift_period_img(vect_of_shift):
    
    image_file = Image.open("/Users/victorstorchan/Desktop/RA/ra/Texture.jpg")
    image_signal = image_file.convert('L')# convert image to black and white
    signal=np.asarray(image_signal)
    signal2=signal[1250:1601,1200:1551]
    (a,b)=signal2.shape
    randommat=np.zeros((a,b),dtype=complex)
    len_shift=len(vect_of_shift)
    signal_noisy=np.zeros((a,b),dtype=complex)
    base=signal2
    #define noisy signal:
    m=signal2.shape[0]
    n=signal2.shape[1]
    for i in range(m):
        for j in range(n):
            randommat[i][j]+=complex(np.random.normal(0,1,1))
    for i in range(m):
        for j in range(n):
            signal_noisy[i][j]= signal2[i][j]#+randommat[i][j]
    
    #signal with shift:
    shifted_signals=[]
    shifted_signals_1=[]
    shifted_signals_2=[]
    #e_1 is vertical, e_2 is horizontal
    max_shift_column=max([b for (a,b) in vect_of_shift])
    max_shift_raw=max([a for (a,b) in vect_of_shift])
    base_period=np.lib.pad(base, ((max_shift_raw, 0),(max_shift_column,0)), 'wrap')
    for s in range(len_shift):
        y=np.zeros((m,n),dtype=complex)
        y1=np.zeros((m,n),dtype=complex)
        y2=np.zeros((m,n),dtype=complex)
        for k in range(m):
            for l in range(n):
                if (k-vect_of_shift[s][0]>=0 and l-vect_of_shift[s][1]>=0 ):
                    y[k][l]=signal_noisy[k-vect_of_shift[s][0]][l-vect_of_shift[s][1]]
                    y1[k][l]=signal_noisy[k-vect_of_shift[s][0]][l-vect_of_shift[s][1]]*np.exp(2J*np.pi*k/m)
                    y2[k][l]=signal_noisy[k-vect_of_shift[s][0]][l-vect_of_shift[s][1]]*np.exp(2J*np.pi*l/n)
                elif (k-vect_of_shift[s][0]<0 and l-vect_of_shift[s][1]>=0):
                    y[k][l]=base_period[max_shift_raw-vect_of_shift[s][0]-1+k][max_shift_column+l-vect_of_shift[s][1]]
                    y1[k][l]=base_period[max_shift_raw-vect_of_shift[s][0]-1+k][max_shift_column+l-vect_of_shift[s][1]]*np.exp(2J*np.pi*k/m)
                    y2[k][l]=base_period[max_shift_raw-vect_of_shift[s][0]-1+k][max_shift_column+l-vect_of_shift[s][1]]*np.exp(2J*np.pi*l/n)
                elif(k-vect_of_shift[s][0]>=0 and l-vect_of_shift[s][1]<0):
                    y[k][l]=base_period[max_shift_raw+k-vect_of_shift[s][0]][max_shift_column-vect_of_shift[s][1]-1+l]
                    y1[k][l]=base_period[max_shift_raw+k-vect_of_shift[s][0]][max_shift_column-vect_of_shift[s][1]-1+l]*np.exp(2J*np.pi*k/m)
                    y2[k][l]=base_period[max_shift_raw+k-vect_of_shift[s][0]][max_shift_column-vect_of_shift[s][1]-1+l]*np.exp(2J*np.pi*l/n)
                else:
                    y[k][l]=base_period[max_shift_raw-vect_of_shift[s][0]-1+k][max_shift_column-vect_of_shift[s][1]-1+l]
                    y1[k][l]=base_period[max_shift_raw-vect_of_shift[s][0]-1+k][max_shift_column-vect_of_shift[s][1]-1+l]*np.exp(2J*np.pi*k/m)
                    y2[k][l]=base_period[max_shift_raw-vect_of_shift[s][0]-1+k][max_shift_column-vect_of_shift[s][1]-1+l]*np.exp(2J*np.pi*l/n)
        shifted_signals.append(y)
        shifted_signals_1.append(y1)
        shifted_signals_2.append(y2)
    print(base_period[max_shift_raw][max_shift_column]==signal2[0][0])    
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
    for s in range(len_shift):
        for i in range(n):
            for k in range(n):
                A_mat1[s][a*i+k]=G1[s][i,k]
                    
    A_mat2=np.zeros((len_shift,n**2),dtype=complex)
    for s in range(len_shift):
        for i in range(n):
            for k in range(n):
                A_mat2[s][a*i+k]=G2[s][i,k]
                    
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
    
    output1=np.zeros(len_shift)
    for i in range(len_shift):
        output1[i]=-n*polar(-v1[i,0])[1]/(2*np.pi).real
    output2=np.zeros(len_shift)
    for i in range(len_shift):
        output2[i]=-m*polar(-v2[i,0])[1]/(2*np.pi).real
    
    #return Image.fromarray(signal2.astype(np.uint8),'L')
    
    return (output1,output2)




#test
shift_period_img([(0,0),(30,100),(20,10),(6,27)])

