#set data

set.seed(1234)
n = 200; p = 100; s = 10
x = matrix(rnorm(n*p),n,p)
beta = c(rnorm(s),rep(0,p-s)) 
y = x %*% beta + rnorm(n)
n0 = 100
x0 = matrix(rnorm(n0*p),n0,p)
y0 = x0 %*% beta + rnorm(n0)
repl=2
funs=lasso.funs()


#test

test_that("Split Error", {
  expect_error(conformal.pred.msplit(x, y, x0, alpha=0.1,
                                     train.fun=funs$train,
                                     predict.fun=funs$predict,
                                     B = repl, split=5), "The 'split' argument should either be NULL or a matrix of dimension B x n1,
         where n1 is the number of observation in the training set !  ")
})

test_that("Correct type", {
  expect_type(conformal.pred.msplit(x, y, x0, alpha=0.1,
                                    train.fun=funs$train,
                                    predict.fun=funs$predict,
                                    B = repl)$lo, "double")
  expect_type(conformal.pred.msplit(x, y, x0, alpha=0.1,
                                    train.fun=funs$train,
                                    predict.fun=funs$predict,
                                    B = repl)$up, "double")
})

test_that("Correct class", {
  expect_is(conformal.pred.msplit(x, y, x0, alpha=0.1,
                                  train.fun=funs$train,
                                  predict.fun=funs$predict,
                                  B = repl)$lo, c("matrix", "array"))
  expect_is(conformal.pred.msplit(x, y, x0, alpha=0.1,
                                  train.fun=funs$train,
                                  predict.fun=funs$predict,
                                  B = repl)$up, c("matrix", "array"))
})

test_that("Correct dimensions", {
  expect_equal(nrow(conformal.pred.msplit(x, y, x0, alpha=0.1,
                                          train.fun=funs$train,
                                          predict.fun=funs$predict,
                                          B = repl)$lo), n0)
  expect_equal(nrow(conformal.pred.msplit(x, y, x0, alpha=0.1,
                                          train.fun=funs$train,
                                          predict.fun=funs$predict,
                                          B = repl)$up), n0)
})

test_that("Silent Code", {
  expect_silent(conformal.pred.msplit(x, y, x0, alpha=0.1,
                                      train.fun=funs$train,
                                      predict.fun=funs$predict,
                                      B = repl))
})

