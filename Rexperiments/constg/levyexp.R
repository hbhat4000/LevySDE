# clear memory
rm(list=ls(all=TRUE))

# load packages
library('fourierin')
library('pracma')

# set parameters
alpha = 1
h = 0.01

f <- function(x)
{
    atan(x)
}

g <- function(x)
{
    1.0
}

p <- function(x)
{
    dnorm(x, mean=0.5, sd=0.1)
}

# set up grids
xmin = -40
xmax = 40
xres = 4096
dx = (xmax-xmin)/xres
xvec = xmin + dx*c(0:(xres-1))

smin = -100
smax = 100
sres = xres
ds = (smax-smin)/sres
svec = smin + ((smax-smin)/sres)*c(0:(sres-1))

bigmats = meshgrid(xvec,svec)
x = bigmats$X
s = bigmats$Y

kernel = exp(1i*s*(x + f(x)*h) - h*(abs(s*g(x))^alpha))

#
# test code to check that meshgrid does what we think it does
#
# testkernel = matrix(0,nrow=sres,ncol=xres)
# for (i in c(1:sres))
# {
#     for (j in c(1:xres))
#     {
#         a = exp(1i*svec[i]*(xvec[j] + f(xvec[j])*h))
#         b = exp(-h*(abs(svec[i]*g(xvec[j])))^alpha)
#         testkernel[i,j] = a*b
#     }
# }

# initialize
den = p(xvec)

# step forward!
numsteps = 100
for (n in c(1:numsteps))
{
    psi = dx*(kernel %*% den)
    print(paste("Normalization: ", psi[sres/2 + 1]))
    psi[sres/2 + 1] = 1
    newden = fourierin(f=psi, lower_int=smin, upper_int=smax, lower_eval=xmin, upper_eval=xmax, const_adj=-1, freq_adj=-1)
    print(paste("L1 norm of Im(newden): ", sum(abs(Im(newden$values)))))
    den = Re(newden$values)
    plot(xvec, den)
}

load('denBR40.RData')
load('nor40.RData')
lines(denBR$x,denBR$y*nor[1]/nor[2],col='red')




