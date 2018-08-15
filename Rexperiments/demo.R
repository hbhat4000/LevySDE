# clear memory
rm(list=ls(all=TRUE))

# load packages
library('fourierin')

######################################################
#                                                    #
# DEMO just so that we remember what fourierin does! #
#                                                    #
######################################################

# define a probability density function
mymean = 0.5
mysd = 0.1
p <- function(x)
{
    dnorm(x, mean=mymean, sd=mysd)
}

# compute characteristic function of the density
testCF = fourierin(f=p, lower_int=-2,upper_int=2,lower_eval=-100,upper_eval=100,resolution=512,const_adj=1,freq_adj=1)

# plot the computed and true characteristic functions
s = testCF$w
trueCF = exp(1i*mymean*s - 0.5*(mysd*s)^2)
par(mfrow=c(2,2))
plot(s, Re(testCF$values))
lines(s, Re(trueCF), col='red')
plot(s, Im(testCF$values))
lines(s, Im(trueCF), col='red')

# now try to recover the original density from the computed characteristic function
testDEN = fourierin(f=testCF$values, lower_int=-100, upper_int=100, lower_eval=-2, upper_eval=2, resolution=512, const_adj=-1, freq_adj=-1)
x = testDEN$w
plot(x, Re(testDEN$values))
lines(x, p(x), col='red')
plot(x, Im(testDEN$values))




