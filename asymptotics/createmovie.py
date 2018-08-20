import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
import scipy.stats as ss

FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='Levy SDE 2', artist='Matplotlib', comment='HSB')
writer = FFMpegWriter(fps=15, metadata=metadata)

# set the value of h
thish = 0.01

# set the initial condition parameters
mymean = 0.5
mysd = 0.1

# set the initial pdf
def initp(x):
    return ss.norm.pdf(x,loc=mymean,scale=mysd)

# set the initial charfun
def initpcf(s):
    return np.exp(1j*mymean*s - 0.5*(mysd*s)**2)

# set up grids of s and u values
smin = -100
smax = 100
sres = 10000
ds = (smax-smin)/sres
svec = smin + ds*np.arange(sres)
uvec = np.copy(svec)

# set up x grid for plotting
Nxplot = 16384
Lplot = 12.8
dxplot = Lplot/Nxplot
xplot = np.arange(-Lplot/2,Lplot/2,dxplot)

# Fourier matrix
fouriermat = np.zeros((Nxplot,sres),dtype='D')
for i in range(Nxplot):
    fouriermat[i,] = np.exp(-1j*uvec*xplot[i])/(2.0*np.pi)

fig = plt.figure()
l, = plt.plot([], [], 'k-.')
plt.xlim(-Lplot/2,Lplot/2)

kernelmat = np.load('kernelmat2.npy')

with writer.saving(fig, "levysde2.mp4", 101):
    pdf = initp(xplot)
    l.set_data(xplot, pdf)
    plt.ylim(-1.0e-4,np.max(pdf)*1.1)
    writer.grab_frame()
    cf = initpcf(svec)
    for i in range(100):
        cf = ds*np.dot(kernelmat, cf)
        pdf = np.real(dxplot*np.dot(fouriermat, cf))
        l.set_data(xplot, pdf)
        plt.ylim(-1.0e-4,np.max(pdf)*1.1)
        writer.grab_frame()



