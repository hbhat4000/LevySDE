# clear memory
rm(list=ls(all=TRUE))

# load packages
library('fourierin')
library('pracma')

# set parameters
alpha = 1
h = 0.01
mymean = 0.5
mysd = 0.1

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
    dnorm(x, mean=mymean, sd=mysd)
}

pcf <- function(s)
{
    exp(1i*mymean*s - 0.5*(mysd*s)^2)
}
    

# set up grids
xmin = -10
xmax = 10
xres = 2048
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

# s is fixed along each row
# x is fixed along each column
kernel = exp(1i*s*(x + f(x)*h) - h*(abs(s*g(x))^alpha))
newkernel = matrix(0,nrow=(sres-1),ncol=sres)
for (i in c(1:(sres-1)))
{
    if (i > sres/2)
        ii = i + 1
    else
        ii = i

    print(ii)
    temp = fourierin(f=kernel[ii,], lower_int=xmin, upper_int=xmax, lower_eval=smin, upper_eval=smax, const_adj=-1, freq_adj=-1)
    newkernel[i,] = temp$values
}

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
psi = pcf(svec)

# step forward!
numsteps = 100
for (n in c(1:numsteps))
{
    psitemp = ds*(newkernel %*% psi)
    psi[1:(sres/2)] = psitemp[1:(sres/2)]
    psi[sres/2 + 1] = 1
    psi[(sres/2 + 2):sres] = psitemp[(sres/2 + 1):(sres-1)]
    den = fourierin(psi, lower_int=smin, upper_int=smax, lower_eval=xmin, upper_eval=xmax, resolution=xres, const_adj=-1, freq_adj=-1)
    plot(xvec, Re(den$values))
}

load('denBR10.RData')
load('nor10.RData')
lines(denBR$x,denBR$y*nor[1]/nor[2],col='red')



