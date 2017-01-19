## Many methods: test out conformal intervals across a variety of
## low-dimensional settings
library(conformalInference)

# Set some overall simulation parameters
n = 200; p = 20 # Numbers of observations and features
s = 5 # Number of truly relevant features
n0 = 100 # Number of points at which to make predictions
nrep = 50 # Number of repetitions for a given setting
sigma = 1 # Marginal error standard deviation
bval = 8 # Magnitude of nonzero coefficients

# Define split conformal inference function for the lasso
my.lasso.funs = lasso.funs(nlambda=50)
my.split.lasso.fun = function(x, y, x0) {
  conformal.pred.split(x,y,x0,alpha=0.1,verb="\t\t",
                      train.fun=my.lasso.funs$train,
                      predict.fun=my.lasso.funs$predict)
}

# Define split conformal inference function for the elastic net 
my.elastic.funs = elastic.funs(gamma=0.8,nlambda=50)
my.split.elastic.fun = function(x, y, x0) {
  conformal.pred.split(x,y,x0,alpha=0.1,verb="\t\t",
                      train.fun=my.elastic.funs$train,
                      predict.fun=my.elastic.funs$predict)
}

# Define split conformal inference function for sparse additive model
my.sam.funs = sam.funs(m=5,nlambda=50,lambda.min.ratio=0.01)
my.split.sam.fun = function(x, y, x0) {
  conformal.pred.split(x,y,x0,alpha=0.1,verb="\t\t",
                      train.fun=my.sam.funs$train,
                      predict.fun=my.sam.funs$predict)
}

# Define split conformal inference function for random forest
my.rf.funs = rf.funs()
my.split.rf.fun = function(x, y, x0) {
  conformal.pred.split(x,y,x0,alpha=0.1,verb="\t\t",
                      train.fun=my.rf.funs$train,
                      predict.fun=my.rf.funs$predict)
}

# Define split conformal inference function for forward stepwise
my.step.funs = lars.funs(max.steps=20)
my.split.step.fun =  function(x, y, x0) {
  conformal.pred.split(x,y,x0,alpha=0.1,verb="\t\t",
                      train.fun=my.step.funs$train,
                      predict.fun=my.step.funs$predict)
}

# Now put together a list with our three conformal inference functions
conformal.pred.funs = vector(mode="list",length=5)
conformal.pred.funs = list(my.split.lasso.fun, my.split.elastic.fun,
  my.split.sam.fun, my.split.rf.fun, my.split.step.fun)
names(conformal.pred.funs) = c("Lasso","Elastic net","Spam",
       "Random forest", "Stepwise")

path = "rds/many.lo."
source("sim.setting.a.R")
source("sim.setting.b.R")
source("sim.setting.c.R")
source("sim.setting.d.R")
#source("sim.setting.e.R")
#source("sim.setting.f.R")

