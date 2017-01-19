## Linear regression: test out conformal intervals and parametric intervals
## across a variety of settings

\dontrun{  
# Set some overall simulation parameters
n = 100; p = 10 # Numbers of observations and features
s = 10 # Number of truly relevant features
n0 = 100 # Number of points at which to make predictions
nrep = 50 # Number of repetitions for a given setting
sigma = 1 # Marginal error standard deviation
alpha = 0.1 # Miscoverage level

# Define conformal inference functions: these are basically just wrappers
# around a particular instatiation of conformal.pred or conformal.pred.split,
# i.e., particular training and prediction functions
my.lm.funs = lm.funs()
my.conf.fun = function(x, y, x0) {
  conformal.pred(x,y,x0,alpha=alpha,verb="\t\t",
                train.fun=my.lm.funs$train,
                predict.fun=my.lm.funs$predict)
}
my.split.fun = function(x, y, x0) {
  conformal.pred.split(x,y,x0,alpha=alpha,
                      train.fun=my.lm.funs$train,
                      predict.fun=my.lm.funs$predict)
}

# Hack together our own "conformal" inference function, really, just one that
# returns the parametric intervals
my.param.fun = function(x, y, x0) {
  lm.obj = lm(y~x)
  mat = predict(lm.obj,list(x=x0),interval="predict",level=1-alpha)
  # Return proper outputs in proper formatting
  return(list(pred=mat[,"fit",drop=F],lo=mat[,"lwr",drop=F],
              up=mat[,"upr",drop=F],fit=matrix(lm.obj$fit,ncol=1)))
}

# Now put together a list with our three conformal inference functions
conformal.pred.funs = list(my.conf.fun, my.split.fun, my.param.fun)
names(conformal.pred.funs) = c("Conformal","Split conformal","Parametric")

# Set some path, where we can save the rds files containing sim results
path = "."

# Simulation setting A: classical setup --- linear (fully supported) mean
# function; normal, homoskedastic errors; and normal, uncorrelated features
cat("========================================\n")
cat("                SETTING A               \n")
cat("========================================\n")
simA = sim.master(n, p, conformal.pred.funs, n0=n0, nrep=nrep, verb=TRUE,
  file=paste0(path,"simA.rds"), x.dist="normal", cor="none",
  mean.fun="linear", s=s, error.dist="normal", sigma=sigma, snr=snr,
  sigma.type="const")

# Simulation setting B: now, nonlinear mean, and t2 errors
cat("========================================\n")
cat("                SETTING B               \n")
cat("========================================\n")
simB = sim.master(n, p, conformal.pred.funs, n0=n0, nrep=nrep, verb=TRUE,
  file=paste0(path,"simB.rds"), x.dist="normal", cor="none",
  mean.fun="additive", m=4, s=s, error.dist="t", df=2, sigma=sigma, snr=snr,
  sigma.type="const")

# Simulation setting C: now, nonlinear mean, heteroskedastic t2 errors, and
# auto-correlated mix of x variables
cat("========================================\n")
cat("                SETTING C               \n")
cat("========================================\n")
simC = sim.master(n, p, conformal.pred.funs, n0=n0, nrep=nrep, verb=TRUE,
  file=paste0(path,"simC.rds"), x.dist="mix", cor="auto", k=5,
  mean.fun="linear", s=s, error.dist="t", df=2, sigma=sigma, snr=snr,
  sigma.type="var")
}
