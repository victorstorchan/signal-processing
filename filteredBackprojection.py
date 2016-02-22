import numpy as np
radon_data= open('/Users/victorstorchan/Desktop/radon_data.dat','r')
radon_data_list=[]
for line in radon_data:
    radon_data_list.append(line.split())
radon_data.close()
radon_data=[]*64
for i in range(len(radon_data_list)):
    radon_data.append(map(float, radon_data_list[i][0].split(',')))
    
"""j'ai mes données dans une matrice radon_data_array, de taille: 64*256"""

radon_data_array=np.asarray(radon_data)



import numpy as np
#from mkfilt import slkernel
import matplotlib.pyplot as plt

# load data
data = np.loadtxt('/Users/victorstorchan/Desktop/radon_data.dat',delimiter=',')
data_poisson= np.ones((64,256))
for i in range(64):
    for j in range(256):
        data_poisson[i][j]=np.random.poisson((4*10**3)*data[i][j],1)
# build t and theta grids
nt     = 256 
ntheta = 64

k = np.linspace(-nt/2,nt/2-1,nt)
t = (k - 0.5)/(255)
dt = 128

l = np.linspace(0,ntheta-1,ntheta)
theta = (np.pi*(l))/(ntheta-1)

# output image grid
nout = 128
m = np.linspace(-nout/2,nout/2-1,nout)
x = (m+0.5)/nout
y = (m+0.5)/nout
dx = x[1] - x[0]
recon = np.zeros([nout,nout])

# building the filter
b=np.pi*dt 
rps=1/b
s = dt*np.linspace(-dt,dt-1,1)
bs = np.linspace(-2*dt,2*dt-1,4*nout)/(dt*rps)
u = np.zeros(len(bs))
u = (np.pi/2 - bs*np.sin(bs))/((np.pi/2)**2 - bs**2)
u = u/(2*np.pi**3)
wb = u/(rps**2)

plt.figure()
plt.plot(wb)

# Plot sinogram
plt.figure()
plt.title('Sinogram')
plt.xlabel('t')
plt.ylabel('theta')
plt.imshow(data)
plt.gray()

# For each angle, loop over all XY's

for ia in range(len(theta)):
    tmp   = np.convolve(data_poisson[ia,:],wb)
    fdata = tmp[2*dt+1:4*dt+1]/dt
    fdata[2*dt-1] = 0.0
    for i in range(len(x)):
        for j in range(len(y)):
            ti  = x[i]*np.cos(theta[ia]) + y[j]*np.sin(theta[ia])
            f   = (ti)*dt
            iff = 2*np.floor(f)
            fx  = f - iff
            gx  = 1. - fx
            idx = int(max(1,iff+dt+1)) # shifting iff so it lies on
            idx = int(min(idx,2*dt-1)-1) # our grid
            #print "x=%f y=%f theta=%f t=%f f=%f it=%f fx=%f gx=%f idx=%d" %(x[i],y[j],theta[ia],ti,f,iff,fx,gx,idx)
            #recon[j][i] = recon[j][i] + data[ia][idx]
            recon[j][i] = recon[j][i] + 0.5*fdata[idx] + 0.5*fdata[idx+1]
        #recon[i][j] = recon[i][j] + gx*data[ia][iff] + fx*data[ia][iff+1] #TODO: improve
#        # with linear interpolation
#
## Plot reconstruction
plt.figure()
plt.title('Filtered Back Projection of poisson distrib.')
plt.xlabel('X')
plt.ylabel('Y')
plt.imshow(recon)
plt.gray()

plt.show()
