# Generate some example training data
set.seed(33)
n = 100; p = 120; s = 10
x = matrix(rnorm(n*p),n,p)
beta = c(rnorm(s),rep(0,p-s)) 
y = x %*% beta + rnorm(n)

# Generate some example test data
n0 = 1
x0 = matrix(rnorm(n0*p),n0,p)
y0 = x0 %*% beta + rnorm(n0)

funs = lm.funs()

# Conformal inference, and jacknife and split conformal versions
out.full.scaled = conformal.quant(x, y, x0, method="scaled", nthreads = 8,
                                      alpha=0.1, 
                                      gamma=0.1/2, verbose=FALSE,w=NULL,
                                      mad.train.fun=NULL,
                                      mad.predict.fun=NULL, num.grid.pts=3, grid.factor=1.25)
  
out.full.classic = conformal.quant(x, y, x0, method="classic", nthreads = 8,
                                               alpha=0.1, 
                                               gamma=0.1/2, verbose=FALSE,w=NULL,
                                               mad.train.fun=NULL,
                                               mad.predict.fun=NULL, num.grid.pts=3, grid.factor=1.25)
  
out.full.median = conformal.quant(x, y, x0, method="median", nthreads = 8,
                                                alpha=0.1, 
                                                gamma=0.1/2, verbose=FALSE,w=NULL,
                                                mad.train.fun=NULL,
                                                mad.predict.fun=NULL, num.grid.pts=3, grid.factor=1.25)



out.split.scaled = conformal.quant.split(x, y, x0, method="scaled", nthreads = 8,
                                                    alpha=0.1, 
                                                    gamma=0.1/2, split=NULL, 
                                                    seed=NULL, verbose=FALSE)
out.split.classic = conformal.quant.split(x, y, x0, method="classic", nthreads = 8,
                                         alpha=0.1, 
                                         gamma=0.1/2, split=NULL, 
                                         seed=NULL, verbose=FALSE)
out.split.median = conformal.quant.split(x, y, x0, method="median", nthreads = 8,
                                         alpha=0.1, 
                                         gamma=0.1/2, split=NULL, 
                                         seed=NULL, verbose=FALSE)
