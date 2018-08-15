# clear memory
rm(list=ls(all=TRUE))

# load required package
library(fourierin)

# define drift and diffusion
myf <- function(x) { atan(x) }
myg <- function(x) { sqrt(1+x^2) }

# time step
h = 0.1

# Levy stable parameter
alpha = 1.5

# phi function
phi <- function(s, y)
{
  exp(1i*s*(y + myf(y)*h))*exp(-h*abs(s*myg(y))^alpha)
}

# set resolution, upper/lower bounds
res = 1024
bnd = 10
mydel = 2*bnd/res
sgrid = seq(from=-bnd,to=bnd,by=mydel)[1:res]
ugrid = sgrid

# build k matrix
kmat = matrix(0, nrow=res, nco=res)
for (i in c(1:res))
{
  thisphi <- function(y) { phi(sgrid[i],y) }
  out <- fourierin(f = thisphi, lower_int = -2*bnd, upper_int = 2*bnd,
                   lower_eval = -bnd, upper_eval = bnd,
                   const_adj = -1, freq_adj = -1,
                   resolution = res)
  kmat[i,] = out$values
}


mysd = 1/2
mypdf <- function(x) { dnorm(x, mean=0, sd=mysd) }
truecf <- function(t) { exp(-0.5*mysd^2*t^2) }
out <- fourierin(f = mypdf, lower_int = -50, upper_int = 50,
                 lower_eval = -bnd, upper_eval = bnd,
                 const_adj = 1, freq_adj = 1,
                 resolution = res)
initcf = out$values





