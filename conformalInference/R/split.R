#' Split conformal prediction intervals.
#'
#' Compute prediction intervals using split conformal inference.
#'
#' @param x Matrix of features, of dimension (say) n x p.
#' @param y Vector of responses, of length (say) n.
#' @param x0 Matrix of features, each row being a point at which we want to
#'   form a prediction interval, of dimension (say) n0 x p. 
#' @param train.fun A function to perform model training, i.e., to produce an
#'   estimator of E(Y|X), the conditional expectation of the response variable
#'   Y given features X. Its input arguments should be x: matrix of features,
#'   and y: vector of responses.
#' @param predict.fun A function to perform prediction for the (mean of the)
#'   responses at new feature values. Its input arguments should be out: output
#'   produced by train.fun, and newx: feature values at which we want to make
#'   predictions.
#' @param alpha Miscoverage level for the prediction intervals, i.e., intervals
#'   with coverage 1-alpha are formed. Default for alpha is 0.1.
#' @param mad.train.fun A function to perform training on the absolute residuals
#'   i.e., to produce an estimator of E(R|X) where R is the absolute residual
#'   R = |Y - m(X)|, and m denotes the estimator produced by train.fun.
#'   This is used to scale the conformal score, to produce a prediction interval
#'   with varying local width. The input arguments to mad.train.fun should be
#'   x: matrix of features, and y: vector of absolute residuals. The default for
#'   mad.train.fun is NULL, which means that no training is done on the absolute
#'   residuals, and the usual (unscaled) conformal score is used. Note that if
#'   mad.train.fun is non-NULL, then so must be mad.predict.fun (see next).
#' @param mad.predict.fun A function to perform prediction for the (mean of the)
#'   absolute residuals at new feature values. Its input arguments should be
#'   out: output produced by mad.train.fun, and newx: feature values at which we
#'   want to make predictions. The default for mad.predict.fun is NULL, which
#'   means that no local scaling is done for the conformal score, i.e., the
#'   usual (unscaled) conformal score is used.
#' @param split Indices that define the data-split to be used (i.e., the indices
#'   define the first half of the data-split, on which the model is trained).
#'   Default is NULL, in which case the split is chosen randomly.
#' @param seed Integer to be passed to set.seed before defining the random 
#'   data-split to be used. Default is NULL, which effectively sets no seed.
#'   If both split and seed are passed, the former takes priority and the latter
#'   is ignored.
#' @param verbose Should intermediate progress be printed out? Default is FALSE.
#'
#' @return A list with the following components: pred, lo, up, fit, split. The
#'   first three are matrices of dimension n0 x m, and the fourth is a matrix of
#'   dimension n x m. Recall that n0 is the number of rows of x0, and m is the
#'   number of tuning parameter values internal to predict.fun. In a sense, each
#'   of the m columns really corresponds to a different prediction function;
#'   see details below. Hence, the rows of the matrices pred, lo, up give
#'   the predicted value, and lower and upper confidence limits (from split
#'   conformal inference), respectively, for the response at the n0 points
#'   given in x0. The rows of fit give the fitted values for the n points given
#'   in x. Finally, split contains the indices used for the first half of the
#'   data-split.
#'
#' @details For concreteness, suppose that we want to use the predictions from
#'   forward stepwise regression at steps 1 through 5 in the path. In this case,
#'   there are m = 5 internal tuning parameter values to predict.fun, in the
#'   notation used above, and each of the returned matrices pred, lo, up, fit
#'   will have 5 columns (one for each step of the forward stepwise path).
#'   The code is structured in this way so that we may defined a single pair of
#'   functions train.fun and predict.fun, over a set of m = 5 tuning parameter
#'   values, instead of calling the conformal function separately m = 5 times.
#'
#' @seealso \code{\link{conformal.pred}}, 
#'   \code{\link{conformal.pred.jack}}, \code{\link{conformal.pred.roo}}
#' @author Ryan Tibshirani
#' @references "Distribution-Free Predictive Inference for Regression" by 
#'   Jing Lei, Max G'Sell, Alessandro Rinaldo, Ryan Tibshirani, and Larry
#'   Wasserman, https://arxiv.org/pdf/1604.04173.pdf, 2016.
#' @example examples/ex.conformal.pred.split.R
#' @export conformal.pred.split

conformal.pred.split = function(x, y, x0, train.fun, predict.fun, alpha=0.1,
  mad.train.fun=NULL, mad.predict.fun=NULL, split=NULL, seed=NULL,
  verbose=FALSE) {

  # Set up data
  x = as.matrix(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)
  x0 = matrix(x0,ncol=p)
  n0 = nrow(x0)
  
  # Check input arguments
  check.args(x=x,y=y,x0=x0,alpha=alpha,train.fun=train.fun,
             predict.fun=predict.fun,mad.train.fun=mad.train.fun,
             mad.predict.fun=mad.predict.fun)

  # Users may pass in a string for the verbose argument
  if (verbose == TRUE) txt = ""
  if (verbose != TRUE && verbose != FALSE) {
    txt = verbose
    verbose = TRUE
  }
  
  # If the user passed indices for the split, use them
  if (!is.null(split)) i1 = split
  # Otherwise make a random split
  else {
    if (!is.null(seed)) set.seed(seed)
    i1 = sample(1:n,floor(n/2))
  }
  i2 = (1:n)[-i1]
  n1 = length(i1)
  n2 = length(i2)
  
  if (verbose) {
    cat(sprintf("%sSplitting data into parts of size %i and %i ...\n",
                txt,n1,n2))
    cat(sprintf("%sTraining on first part ...\n",txt))
  }
  
  # Train on first part
  out = train.fun(x[i1,,drop=F],y[i1])
  fit = matrix(predict.fun(out,x),nrow=n)
  pred = matrix(predict.fun(out,x0),nrow=n0)
  m = ncol(pred)
  
  if (verbose) {
    cat(sprintf("%sComputing residuals and quantiles on second part ...\n",txt))
  }
  
  # Get residuals and quantiles on second
  res = abs(y[i2] - matrix(predict.fun(out,x[i2,,drop=F]),nrow=n2))
  lo = up = matrix(0,n0,m)
  
  for (l in 1:m) {
    # Local scoring?
    if (!is.null(mad.train.fun) && !is.null(mad.predict.fun)) {
      res.train = abs(y[i1] - fit[i1,l])
      mad.out = mad.train.fun(x[i1,,drop=F],res.train)
      mad.x2 = mad.predict.fun(mad.out,x[i2,,drop=F])
      mad.x0 = mad.predict.fun(mad.out,x0)
      res[,l] = res[,l] / mad.x2
    }
    else {
      mad.x0 = rep(1,n0)
    }
    
    q = conformal.quantile(sort(res[,l]),alpha)
    lo[,l] = pred[,l] - q * mad.x0
    up[,l] = pred[,l] + q * mad.x0
  }
 
  return(list(pred=pred,lo=lo,up=up,fit=fit,split=i1))
}
