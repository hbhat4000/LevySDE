library(fourierin)

mysd = 1/20

mypdf <- function(x) { dnorm(x, mean=0, sd=mysd) }

truecf <- function(t) { exp(-0.5*mysd^2*t^2) }

out <- fourierin(f = mypdf, lower_int = -5, upper_int = 5,
                 lower_eval = -100, upper_eval = 100,
                 const_adj = 1, freq_adj = 1,
                 resolution = 800)

plot(out$w,Re(out$values))
lines(out$w,truecf(out$w),col='red')
