#set data

set.seed(1234)
n = 200; p = 100; s = 10
x = matrix(rnorm(n*p),n,p)
beta = c(rnorm(s),rep(0,p-s)) 
y = x %*% beta + rnorm(n)
n0 = 2
x0 = matrix(rnorm(n0*p),n0,p)
y0 = x0 %*% beta + rnorm(n0)
funs=lasso.funs()


#test

test_that("Argument Errors", {
  expect_error(conformal.quant(x, y, x0, method="yax", nthreads = 8,
                                         alpha=0.1, 
                                         gamma=0.1/2, verbose=FALSE,w=NULL,
                                         mad.train.fun=NULL,
                                         mad.predict.fun=NULL, num.grid.pts=2, grid.factor=1.25), "The 'method' argument is not correct. \n             Please select one of the following:classic, median, scaled.")
  
  expect_error(conformal.quant(x="tfvyr", y, x0, method="classic", nthreads = 8,
                               alpha=0.1, 
                               gamma=0.1/2, verbose=FALSE,w=NULL,
                               mad.train.fun=NULL,
                               mad.predict.fun=NULL, num.grid.pts=2, grid.factor=1.25), "x must be a numeric matrix")
  
  
  expect_error(conformal.quant(x, y, x0="tfvy3r", method="classic", nthreads = 8,
                               alpha=0.1, 
                               gamma=0.1/2, verbose=FALSE,w=NULL,
                               mad.train.fun=NULL,
                               mad.predict.fun=NULL, num.grid.pts=2, grid.factor=1.25), "x0 must be a numeric matrix")
  
})

test_that("Correct type", {
  expect_type(conformal.quant(x, y, x0, method="classic", nthreads = 8,
                              alpha=0.1, 
                              gamma=0.1/2, verbose=FALSE,w=NULL,
                              mad.train.fun=NULL,
                              mad.predict.fun=NULL, num.grid.pts=100, grid.factor=1.25)$lo, "double")
  expect_type(conformal.quant(x, y, x0, method="classic", nthreads = 8,
                              alpha=0.1, 
                              gamma=0.1/2, verbose=FALSE,w=NULL,
                              mad.train.fun=NULL,
                              mad.predict.fun=NULL, num.grid.pts=100, grid.factor=1.25)$up, "double")
  expect_type(conformal.quant(x, y, x0, method="median", nthreads = 8,
                              alpha=0.1, 
                              gamma=0.1/2, verbose=FALSE,w=NULL,
                              mad.train.fun=NULL,
                              mad.predict.fun=NULL, num.grid.pts=100, grid.factor=1.25)$lo, "double")
  expect_type(conformal.quant(x, y, x0, method="median", nthreads = 8,
                              alpha=0.1, 
                              gamma=0.1/2, verbose=FALSE,w=NULL,
                              mad.train.fun=NULL,
                              mad.predict.fun=NULL, num.grid.pts=100, grid.factor=1.25)$up, "double")
  expect_type(conformal.quant(x, y, x0, method="scaled", nthreads = 8,
                              alpha=0.1, 
                              gamma=0.1/2, verbose=FALSE,w=NULL,
                              mad.train.fun=NULL,
                              mad.predict.fun=NULL, num.grid.pts=100, grid.factor=1.25)$lo, "double")
  expect_type(conformal.quant(x, y, x0, method="scaled", nthreads = 8,
                              alpha=0.1, 
                              gamma=0.1/2, verbose=FALSE,w=NULL,
                              mad.train.fun=NULL,
                              mad.predict.fun=NULL, num.grid.pts=100, grid.factor=1.25)$up, "double")
})

test_that("Correct length", {
  expect_length(conformal.quant.split(x, y, x0, method="classic", nthreads = 8,
                                      alpha=0.1, 
                                      gamma=0.1/2, split=NULL, 
                                      seed=NULL, verbose=FALSE)$lo,n0)
  expect_length(conformal.quant.split(x, y, x0, method="classic", nthreads = 8,
                                      alpha=0.1, 
                                      gamma=0.1/2, split=NULL, 
                                      seed=NULL, verbose=FALSE)$up,n0)
  
  expect_length(conformal.quant.split(x, y, x0, method="median", nthreads = 8,
                                      alpha=0.1, 
                                      gamma=0.1/2, split=NULL, 
                                      seed=NULL, verbose=FALSE)$lo,n0)
  expect_length(conformal.quant.split(x, y, x0, method="median", nthreads = 8,
                                      alpha=0.1, 
                                      gamma=0.1/2, split=NULL, 
                                      seed=NULL, verbose=FALSE)$up,n0)
  
  expect_length(conformal.quant.split(x, y, x0, method="scaled", nthreads = 8,
                                      alpha=0.1, 
                                      gamma=0.1/2, split=NULL, 
                                      seed=NULL, verbose=FALSE)$lo,n0)
  expect_length(conformal.quant.split(x, y, x0, method="scaled", nthreads = 8,
                                      alpha=0.1, 
                                      gamma=0.1/2, split=NULL, 
                                      seed=NULL, verbose=FALSE)$up ,n0)
})





