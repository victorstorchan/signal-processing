"""example of how the svd method of python is shifting all the value with respect to the first entry:
During the svd decomposition, w and w* are splitted between  v sigma and v* where v[0]=1 (w=v/c and w*=v*c)""""


import numpy as np
#l1=3, l2=0 l3=5
r=np.matrix([[ 1,np.exp(-2J*np.pi*(3-0)/10),np.exp(-2J*np.pi*(3-5)/10)],
        [np.exp(-2J*np.pi*(0-3)/10),1,np.exp(-2J*np.pi*(0-5)/10)],
        [np.exp(-2J*np.pi*(5-3)/10),np.exp(-2J*np.pi*(5-0)/10),1]])


(V,sigma,V_star)=np.linalg.svd(r,full_matrices=0)
v1=V_star[0].getH()
output=np.zeros(3,dtype=complex)
for i in range(3):
    output[i]=-10*polar(-v1[i,0])[1]/(2*np.pi)
output

"""the result is array([-0.+0.j, -3.+0.j,  2.+0.j]) so when the first shift is 
not 0, the output vector of shifts is shifted by the first shift. This is just
due to the implementation of python's SVD like this example is showing it"""
