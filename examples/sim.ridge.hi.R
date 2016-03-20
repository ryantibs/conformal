## Linear regression: test out conformal intervals and parametric intervals
## across a variety of high-dimensional settings
library(conformalInference)

# Set some overall simulation parameters
n = 500; p = 490 # Numbers of observations and features
s = 10 # Number of truly relevant features
n0 = 100 # Number of points at which to make predictions
nrep = 50 # Number of repetitions for a given setting
sigma = 1 # Marginal error standard deviation
snr = 1 # Signal-to-noise ratio

# Define conformal inference functions: these are basically just wrappers
# around a particular instatiation of conformal.pred or conformal.pred.split,
# i.e., particular training and prediction functions
my.lm.funs = lm.funs()
my.conf.fun = function(x, y, x0) {
  conformal.pred(x,y,x0,alpha=0.1,verb="\t\t",
                train.fun=my.lm.funs$train,
                predict.fun=my.lm.funs$predict)
}
my.split.fun = function(x, y, x0) {
  conformal.pred.split(x,y,x0,alpha=0.1,
                      train.fun=my.lm.funs$train,
                      predict.fun=my.lm.funs$predict)
}

# Hack together our own "conformal" inference function, really, just one that
# returns the parametric intervals
my.param.fun = function(x, y, x0) {
  lm.obj = lm(y~x)
  mat = predict(lm.obj,list(x=x0),interval="predict",level=0.9)
  # Return proper outputs in proper formatting
  return(list(pred=mat[,"fit",drop=F],lo=mat[,"lwr",drop=F],
              up=mat[,"upr",drop=F],fit=matrix(lm.obj$fit,ncol=1)))
}

# Now put together a list with our three conformal inference functions
conformal.pred.funs = list(my.conf.fun, my.split.fun, my.param.fun)
names(conformal.pred.funs) = c("Conformal","Split conformal","Parametric")

conformal.pred.funs = conformal.pred.funs[c(1,3)]
path = "rds/lm.hi."
source("sim.setting.a.R")
source("sim.setting.b.R")
source("sim.setting.c.R")

