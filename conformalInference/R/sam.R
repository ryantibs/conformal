#' Sparse additive models training and prediction functions.
#'
#' Construct training and prediction functions for sparse additive models,
#'   based on the \code{\link{SAM}} package, over a sequence of (given or
#'   internally computed) lambda values.
#'
#' @param m Number of spline functions to be used in the additive expansion of 
#'   each variable. Default is 3.
#' @param lambda Sequence of lambda values over which training is performed.
#'   This must be in decreasing order, and --- this argument should be used with
#'   caution! When used, it is usually best to grab the sequence constructed by
#'   one intial call to SAM (see examples). Default is NULL, which means that
#'   the nlambda, lambda.min.ratio arguments will define the lambda sequence
#'   (see next).
#' @param nlambda Number of lambda values over which training is performed. In
#'   particular, the lambda sequence is defined by nlambda log-spaced values
#'   between lambda.max and lambda.min.ratio * lambda.max, where lambda.max is
#'   the smallest value of lambda at which the solution has all zero components,
#'   and lambda.min.ratio is a small fraction (see next). Default is 30.
#'
#' @return A list with three components: train.fun, predict.fun, active.fun.
#'   The third function is designed to take the output of train.fun, and
#'   reports which features are active for each fitted model contained in
#'   this output. 
#'
#' @details This function is based on the packages \code{\link{SAM}} and
#'   \code{\link{plyr}}. If these packages are not installed, then the function
#'   will abort. 
#'
#' @author Ryan Tibshirani
#' @example examples/ex.sam.funs.R
#' @export sam.funs

sam.funs = function(m=3, lambda=NULL, nlambda=30, lambda.min.ratio=5e-3) {
  # Check for SAM
  if (!require("SAM",quietly=TRUE)) {
    stop("Package SAM not installed (required here)!")
  }
  # Check for plyr
  if (!require("plyr",quietly=TRUE)) {
    stop("Package plyr not installed (required here)!")
  }

  # Check arguments
  check.pos.int(m)
  if (is.null(lambda) && (length(nlambda)!= 1 ||
                          round(nlambda) != nlambda ||
                          nlambda < 1 || nlambda > 100)) {
    stop("nlambda must be an integer between 1 and 100")
  }
  if (!is.null(lambda) && (!is.numeric(lambda) || min(lambda)<=0 ||
                           (order(lambda) != length(lambda):1))) {
    stop("lambda must be a decreasing sequence of positive numbers")
  }
  check.pos.num(lambda.min.ratio)

  train.fun = function(x,y,out=NULL) {
    return(samQL(x,y,p=m,lambda=lambda,nlambda=nlambda,
                 lambda.min.ratio=lambda.min.ratio))
  }

  predict.fun = function(out,newx) {
    return(predict(out,newx)$values)
  }

  active.fun = function(out) {
    return(alply(out$func,2,.fun=function(v) which(v!=0)))
  }
  
  return(list(train.fun=train.fun, predict.fun=predict.fun,
              active.fun=active.fun))
}
