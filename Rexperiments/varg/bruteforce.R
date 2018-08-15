# clear memory
rm(list=ls(all=TRUE))

# set parameters
alpha = 1
h = 1e-2
numsteps = 1e2
mymean = 0.5
mysd = 0.1

f <- function(x)
{
    atan(x)
}

g <- function(x)
{
    sqrt(1 + x^2)
}

numpaths = 1e7

# euler-maruyama
require(foreach)
require(doParallel)
ncores = 24
c1 = makeCluster(ncores)
registerDoParallel(c1)

createsamples <- function()
{
    # initialize
    x = rnorm(n=numpaths, mean=mymean, sd=mysd)

    # step forward
    for (n in c(1:numsteps))
    {
        x = x + f(x)*h + g(x)*rcauchy(n=numpaths,location=0,scale=h)
    }
    x
}

x = foreach(i=1:ncores,.combine='c') %dopar% createsamples()
save(x,file='x5.RData')


