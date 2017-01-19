## Linear regression: test out conformal intervals and parametric intervals
## across a variety of low-dimensional settings
library(conformalInference)

# Set some overall simulation parameters
n = 100; p = 10 # Numbers of observations and features
s = 10 # Number of truly relevant features
n0 = 100 # Number of points at which to make predictions
nrep = 50 # Number of repetitions for a given setting
sigma = 1 # Marginal error standard deviation
bval = 1 # Magnitude of nonzero coefficients
lambda = 0 # Lambda values to try in ridge regression
alpha = 0.1 # Miscoverage level

# Define conformal inference functions: these are basically just wrappers
# around a particular instatiation of conformal.pred, conformal.pred.jack, or
# conformal.pred.split
my.lm.funs = lm.funs(lambda=lambda)
my.conf.fun = function(x, y, x0) {
  conformal.pred(x,y,x0,alpha=alpha,verb="\t\t",
                train.fun=my.lm.funs$train,
                predict.fun=my.lm.funs$predict)
}
my.jack.fun = function(x, y, x0) {
  conformal.pred.jack(x,y,x0,alpha=alpha,verb="\t\t",
                     train.fun=my.lm.funs$train,
                     predict.fun=my.lm.funs$predict,
                     special.fun=my.lm.funs$special)
}
my.split.fun = function(x, y, x0) {
  conformal.pred.split(x,y,x0,alpha=alpha,
                      train.fun=my.lm.funs$train,
                      predict.fun=my.lm.funs$predict)
}

# Hack together our own "conformal" inference function, really, just one that
# returns the parametric intervals
my.param.fun = function(x, y, x0) {
  n = nrow(x); n0 = nrow(x0)
  out = my.lm.funs$train(x,y)
  fit = matrix(my.lm.funs$predict(out,x),nrow=n)
  pred = matrix(my.lm.funs$predict(out,x0),nrow=n0)
  m = ncol(pred)
  
  x1 = cbind(rep(1,n0),x0)
  z = qnorm(1-alpha/2)
  lo = up = matrix(0,n0,m)
  
  for (j in 1:m) {
    sig.hat = sqrt(sum((y - fit[,j])^2)/(n-ncol(x1)))
    g = diag(x1 %*% chol.solve(out$chol.R[[j]], t(x1)))
    lo[,j] = pred[,j] - sqrt(1+g)*sig.hat*z
    up[,j] = pred[,j] + sqrt(1+g)*sig.hat*z
  }
  
  # Return proper outputs in proper formatting
  return(list(pred=pred,lo=lo,up=up,fit=fit))
}

# Now put together a list with all of our conformal inference functions
conformal.pred.funs = list(my.conf.fun, my.jack.fun, my.split.fun, my.param.fun)
names(conformal.pred.funs) = c("Conformal","Jackknife","Split conformal",
       "Parametric")

path = "rds/lm.lo."
source("sim.setting.a.R")
source("sim.setting.b.R")
source("sim.setting.c.R")
