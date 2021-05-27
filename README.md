# KRR-Derivatives
Using Plug-in Kernel ridge regression estimators to estimate function derivatives 

```R
library(Rcpp)
library(RandomFieldsUtils)
library(latex2exp)
source("SourceCode.R")

# True regression function
f0 = function(x){
  #exp(-4*(1-2*x)^2)*(1-2*x)
  sin(8*x)+cos(8*x)+log(4/3+x)
}

# Derivative of true regression function
f0_prime = function(x){
  #-2*exp(-4*(1-2*x)^2)+16*exp(-4*(1-2*x)^2)*(1-2*x)^2
  8*cos(8*x)-8*sin(8*x)+1/(4/3+x)
}

set.seed(1)
n = 500 # number of data points
x = sort(runif(n, 0, 1)) # random design
I = 1:9999
y0 = sapply(x, f0) # true value at samples
# sigma0 = sqrt(0.1^2)
sigma0 = sqrt(0.2^2) # true regression standard deviation
y = y0 + sigma0 * rnorm(n)  #generate some data

n_new = 100 # number of grid points
x_new = seq(0, 1,, n_new) # grid points
f0_new = sapply(x_new, f0) # true value at grid points
f0_new_prime = sapply(x_new, f0_prime) # true derivative at grid points

plot(x_new, f0_new_prime, type = "l", lwd = 5, xlab = "x", ylab = TeX("f'_{02}(x),\\hat{f}'_{02}(x)"), ylim = c(-15,15))

get_GPR("Sobolev")
get_GPR("Matern")

#install.pachakes("locpol")
library(locpol)
bw = regCVBwSelC(x, y, 2, kernel=gaussK, interval = c(0,1))
lpest2 = locPolSmootherC(x, y, x_new, bw, 2, gaussK)
lines(x_new, lpest2$beta1, lty = 2, lwd = 5, col = "blue")

#install.pachakes("pspline")
library(pspline)
ps_y = as.vector(predict(sm.spline(x, y, norder = 2), x_new, 1))
lines(x_new, ps_y, lty = 2, lwd = 5, col = "yellow")
```

