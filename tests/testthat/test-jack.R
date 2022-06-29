#set data

set.seed(1234)
n = 10; p = 12; s = 10
x = matrix(rnorm(n*p),n,p)
beta = c(rnorm(s),rep(0,p-s)) 
y = x %*% beta + rnorm(n)
n0 = 2
x0 = matrix(rnorm(n0*p),n0,p)
y0 = x0 %*% beta + rnorm(n0)
funs=lasso.funs()


#test

test_that("Split Error", {
  expect_error(conformal.pred.jack(x, y, x0, train.fun=funs$train,
                                   predict.fun=funs$predict, alpha=0.1, 
                                              special.fun=NULL, mad.train.fun=NULL, 
                                   mad.predict.fun=NULL, 
                                              verbose=FALSE, plus=5),
               "The argument 'plus' should be a logical \n")
})


test_that("Correct type", {
  expect_type(conformal.pred.jack(x, y, x0, train.fun=funs$train,
                                  predict.fun=funs$predict, alpha=0.1, 
                                  special.fun=NULL, mad.train.fun=NULL, 
                                  mad.predict.fun=NULL, 
                                  verbose=FALSE, plus=FALSE)$lo, "double")
  expect_type(conformal.pred.jack(x, y, x0, train.fun=funs$train,
                                  predict.fun=funs$predict, alpha=0.1, 
                                  special.fun=NULL, mad.train.fun=NULL, 
                                  mad.predict.fun=NULL, 
                                  verbose=FALSE, plus=FALSE)$up, "double")
  expect_type(conformal.pred.jack(x, y, x0, train.fun=funs$train,
                                  predict.fun=funs$predict, alpha=0.1, 
                                  special.fun=NULL, mad.train.fun=NULL, 
                                  mad.predict.fun=NULL, 
                                  verbose=FALSE, plus=FALSE)$fit, "double")
  expect_type(conformal.pred.jack(x, y, x0, train.fun=funs$train,
                                  predict.fun=funs$predict, alpha=0.1, 
                                  special.fun=NULL, mad.train.fun=NULL, 
                                  mad.predict.fun=NULL, 
                                  verbose=FALSE, plus=FALSE)$pred, "double")
  expect_type(conformal.pred.jack(x, y, x0, train.fun=funs$train,
                                  predict.fun=funs$predict, alpha=0.1, 
                                  special.fun=NULL, mad.train.fun=NULL, 
                                  mad.predict.fun=NULL, 
                                  verbose=FALSE, plus=TRUE)$lo, "double")
  expect_type(conformal.pred.jack(x, y, x0, train.fun=funs$train,
                                  predict.fun=funs$predict, alpha=0.1, 
                                  special.fun=NULL, mad.train.fun=NULL, 
                                  mad.predict.fun=NULL, 
                                  verbose=FALSE, plus=TRUE)$up, "double")
  expect_type(conformal.pred.jack(x, y, x0, train.fun=funs$train,
                                  predict.fun=funs$predict, alpha=0.1, 
                                  special.fun=NULL, mad.train.fun=NULL, 
                                  mad.predict.fun=NULL, 
                                  verbose=FALSE, plus=TRUE)$fit, "double")
  expect_type(conformal.pred.jack(x, y, x0, train.fun=funs$train,
                                  predict.fun=funs$predict, alpha=0.1, 
                                  special.fun=NULL, mad.train.fun=NULL, 
                                  mad.predict.fun=NULL, 
                                  verbose=FALSE, plus=TRUE)$pred, "double")
})

test_that("Correct class", {
  expect_is(conformal.pred.jack(x, y, x0, train.fun=funs$train,
                                predict.fun=funs$predict, alpha=0.1, 
                                special.fun=NULL, mad.train.fun=NULL, 
                                mad.predict.fun=NULL, 
                                verbose=FALSE, plus=FALSE)$lo, c("matrix", "array"))
  expect_is(conformal.pred.jack(x, y, x0, train.fun=funs$train,
                                predict.fun=funs$predict, alpha=0.1, 
                                special.fun=NULL, mad.train.fun=NULL, 
                                mad.predict.fun=NULL, 
                                verbose=FALSE, plus=FALSE)$up, c("matrix", "array"))
  expect_is(conformal.pred.jack(x, y, x0, train.fun=funs$train,
                                predict.fun=funs$predict, alpha=0.1, 
                                special.fun=NULL, mad.train.fun=NULL, 
                                mad.predict.fun=NULL, 
                                verbose=FALSE, plus=FALSE)$pred, c("matrix", "array"))
  expect_is(conformal.pred.jack(x, y, x0, train.fun=funs$train,
                                predict.fun=funs$predict, alpha=0.1, 
                                special.fun=NULL, mad.train.fun=NULL, 
                                mad.predict.fun=NULL, 
                                verbose=FALSE, plus=FALSE)$fit, c("matrix", "array"))
  
  expect_is(conformal.pred.jack(x, y, x0, train.fun=funs$train,
                                predict.fun=funs$predict, alpha=0.1, 
                                special.fun=NULL, mad.train.fun=NULL, 
                                mad.predict.fun=NULL, 
                                verbose=FALSE, plus=TRUE)$lo, c("matrix", "array"))
  expect_is(conformal.pred.jack(x, y, x0, train.fun=funs$train,
                                predict.fun=funs$predict, alpha=0.1, 
                                special.fun=NULL, mad.train.fun=NULL, 
                                mad.predict.fun=NULL, 
                                verbose=FALSE, plus=TRUE)$up, c("matrix", "array"))
  expect_is(conformal.pred.jack(x, y, x0, train.fun=funs$train,
                                predict.fun=funs$predict, alpha=0.1, 
                                special.fun=NULL, mad.train.fun=NULL, 
                                mad.predict.fun=NULL, 
                                verbose=FALSE, plus=TRUE)$pred, c("matrix", "array"))
  expect_is(conformal.pred.jack(x, y, x0, train.fun=funs$train,
                                predict.fun=funs$predict, alpha=0.1, 
                                special.fun=NULL, mad.train.fun=NULL, 
                                mad.predict.fun=NULL, 
                                verbose=FALSE, plus=TRUE)$fit, c("matrix", "array"))
})

test_that("Correct dimensions", {
  expect_equal(nrow(conformal.pred.jack(x, y, x0, train.fun=funs$train,
                                        predict.fun=funs$predict, alpha=0.1, 
                                        special.fun=NULL, mad.train.fun=NULL, 
                                        mad.predict.fun=NULL, 
                                        verbose=FALSE, plus=FALSE)$lo), n0)
  expect_equal(nrow(conformal.pred.jack(x, y, x0, train.fun=funs$train,
                                        predict.fun=funs$predict, alpha=0.1, 
                                        special.fun=NULL, mad.train.fun=NULL, 
                                        mad.predict.fun=NULL, 
                                        verbose=FALSE, plus=FALSE)$up), n0)
  expect_equal(nrow(conformal.pred.jack(x, y, x0, train.fun=funs$train,
                                        predict.fun=funs$predict, alpha=0.1, 
                                        special.fun=NULL, mad.train.fun=NULL, 
                                        mad.predict.fun=NULL, 
                                        verbose=FALSE, plus=FALSE)$pred), n0)
  expect_equal(nrow(conformal.pred.jack(x, y, x0, train.fun=funs$train,
                                        predict.fun=funs$predict, alpha=0.1, 
                                        special.fun=NULL, mad.train.fun=NULL, 
                                        mad.predict.fun=NULL, 
                                        verbose=FALSE, plus=FALSE)$fit), n)
  
  
  expect_equal(nrow(conformal.pred.jack(x, y, x0, train.fun=funs$train,
                                        predict.fun=funs$predict, alpha=0.1, 
                                        special.fun=NULL, mad.train.fun=NULL, 
                                        mad.predict.fun=NULL, 
                                        verbose=FALSE, plus=TRUE)$lo), n0)
  expect_equal(nrow(conformal.pred.jack(x, y, x0, train.fun=funs$train,
                                        predict.fun=funs$predict, alpha=0.1, 
                                        special.fun=NULL, mad.train.fun=NULL, 
                                        mad.predict.fun=NULL, 
                                        verbose=FALSE, plus=TRUE)$up), n0)
  expect_equal(nrow(conformal.pred.jack(x, y, x0, train.fun=funs$train,
                                        predict.fun=funs$predict, alpha=0.1, 
                                        special.fun=NULL, mad.train.fun=NULL, 
                                        mad.predict.fun=NULL, 
                                        verbose=FALSE, plus=TRUE)$pred), n0)
  expect_equal(nrow(conformal.pred.jack(x, y, x0, train.fun=funs$train,
                                        predict.fun=funs$predict, alpha=0.1, 
                                        special.fun=NULL, mad.train.fun=NULL, 
                                        mad.predict.fun=NULL, 
                                        verbose=FALSE, plus=TRUE)$fit), n)
})

test_that("Silent Code", {
  expect_silent(conformal.pred.jack(x, y, x0, train.fun=funs$train,
                                    predict.fun=funs$predict, alpha=0.1, 
                                    special.fun=NULL, mad.train.fun=NULL, 
                                    mad.predict.fun=NULL, 
                                    verbose=FALSE, plus=TRUE))
  expect_silent(conformal.pred.jack(x, y, x0, train.fun=funs$train,
                                    predict.fun=funs$predict, alpha=0.1, 
                                    special.fun=NULL, mad.train.fun=NULL, 
                                    mad.predict.fun=NULL, 
                                    verbose=FALSE, plus=FALSE))
})

test_that("Verbose Code", {
  expect_output(conformal.pred.jack(x, y, x0, train.fun=funs$train,
                                    predict.fun=funs$predict, alpha=0.1, 
                                    special.fun=NULL, mad.train.fun=NULL, 
                                    mad.predict.fun=NULL, 
                                    verbose=TRUE, plus=TRUE),
                "%sInitial training on full data set ...\n")
  expect_output(conformal.pred.jack(x, y, x0, train.fun=funs$train,
                                    predict.fun=funs$predict, alpha=0.1, 
                                    special.fun=NULL, mad.train.fun=NULL, 
                                    mad.predict.fun=NULL, 
                                    verbose=TRUE, plus=FALSE),
                "%sInitial training on full data set ...\n")
})
