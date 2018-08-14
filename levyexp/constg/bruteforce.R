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
    1.0
}

numpaths = 1e7

# initialize
x = rnorm(n=numpaths, mean=mymean, sd=mysd)

# euler-maruyama
require(foreach)
require(doParallel)
c1 = makeCluster(24)
registerDoParallel(c1)

ncores = 24
for (n in c(1:numsteps))
{
    print(n)
    myrands = foreach(i=1:24,.combine='c') %dopar% rcauchy(n=numpaths,location=0,scale=h)
    x = x + f(x)*h + g(x)*myrands 
}

plot(density(x))


