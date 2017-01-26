#' Jackknife conformal prediction intervals.
#'
#' Compute prediction intervals using a jackknife variant of conformal
#'   inference.
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
#' @param special.fun A function to compute leave-one-out fitted values in an
#'   an efficient manner. For example, for many linear smoothers (e.g., linear
#'   regression smoothing splines, kernel smoothing --- essentially, anything
#'   for which the fast leave-one-out cross-validation formula holds), we have
#'   the "magic" property:
#'   \deqn{ y_i - \hat{f}^{(-i)}(x_i) = \frac{y_i - \hat{f}(x_i)}{1-S_{ii}}, }
#'   where \eqn{\hat{f}} is the linear smoother trained on the entire data set
#'   \eqn{\hat{f}^{(-i)}} is the linear smoother trained on all but the ith data
#'   point, and \eqn{S_{ii}} is the ith diagonal element of the smoothing
#'   matrix. The input arguments to special.fun should be x: matrix of features, 
#'   y: vector of responses, and out: output produced by train.fun when run on
#'   x, y. The default for special.fun is NULL, which means that no special
#'   leave-one-out formula is leverage, and the usual (manual) way of computing
#'   leave-one-out fitted values is used. Note that if mad.train.fun and
#'   mad.predict.fun are passed, in order to obtain prediction intervals with
#'   locally varying width, then special.fun will be ignored (since the fast
#'   shortcut it offers can no longer be used).
#' @param mad.train.fun A function to perform training on the absolute residuals
#'   i.e., to produce an estimator of E(R|X) where R is the absolute residual
#'   R = |Y - m(X)|, and m denotes the estimator produced by train.fun.
#'   This is used to scale the conformal score, to produce a prediction interval
#'   with varying local width. The input arguments to mad.train.fun should be
#'   x: matrix of features, and y: vector of absolute residuals. The default for
#'   mad.train.fun is NULL, which means that no training is done on the absolute
#'   residuals, and the usual (unscaled) conformal score is used. Note that if
#'   mad.train.fun is non-NULL, then so must be mad.predict.fun (next).
#' @param mad.predict.fun A function to perform prediction for the (mean of the)
#'   absolute residuals at new feature values. Its input arguments should be
#'   out: output produced by mad.train.fun, and newx: feature values at which we
#'   want to make predictions. The default for mad.predict.fun is NULL, which
#'   means that no local scaling is done for the conformal score, i.e., the
#'   usual (unscaled) conformal score is used.
#' @param alpha Miscoverage level for the prediction intervals, i.e., intervals
#'   with coverage 1-alpha are formed. Default for alpha is 0.1.
#' @param verbose Should intermediate progress be printed out? Default is FALSE.
#'
#' @return A list with the following components: pred, lo, up, fit. The first
#'   three are matrices of dimension m x n0, and the last is a matrix of
#'   dimension m x n. Recall that n0 is the number of rows of x0, and m is the
#'   number of tuning parameter values internal to predict.fun. In a sense, each
#'   of the m columns really corresponds to a different prediction function;
#'   see details below. Hence, the rows of the matrices pred, lo, up give
#'   the predicted value, and lower and upper confidence limits, respectively,
#'   for the response at the n0 points given in x0. The predicted values
#'   are computed based on the entire training set, and the confidence limits
#'   are computed using jackknife conformal inference. The rows of fit give the
#'   fitted values for the n points given in x.
#'
#' @details For concreteness, suppose that we want to use the predictions from
#'   forward stepwise regression at steps 1 through 5 in the path. In this case,
#'   there are m = 5 internal tuning parameter values to predict.fun, in the
#'   notation used above, and each of the returned matrices pred, lo, and up will
#'   have 5 rows (one corresponding to each step of the forward stepwise path).
#'   The code is structured in this way so that we may defined a single pair of
#'   functions train.fun and predict.fun, over a set of m = 5 tuning parameter
#'   values, instead of calling the conformal function separately m = 5 times.
#'
#' @seealso \code{\link{conformal.pred}}, 
#'   \code{\link{conformal.pred.split}}, \code{\link{conformal.pred.roo}}
#' @author Ryan Tibshirani
#' @references "Distribution-Free Predictive Inference for Regression" by 
#'   Jing Lei, Max G'Sell, Alessandro Rinaldo, Ryan Tibshirani, and Larry
#'   Wasserman, https://arxiv.org/pdf/1604.04173.pdf, 2016.
#' @references \url{http://www.stat.cmu.edu}
#' @examples ## See examples for conformal.pred function
#' @export conformal.pred.jack

conformal.pred.jack = function(x, y, x0, train.fun, predict.fun, alpha=0.1, 
  special.fun=NULL, mad.train.fun=NULL, mad.predict.fun=NULL, verbose=FALSE) {

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
             mad.predict.fun=mad.predict.fun,special.fun=special.fun)
  
  # Users may pass in a string for the verbose argument
  if (verbose == TRUE) txt = ""
  if (verbose != TRUE && verbose != FALSE) {
    txt = verbose
    verbose = TRUE
  }

  if (verbose) cat(sprintf("%sInitial training on full data set ...\n",txt))
  
  # Train, fit, and predict on full data set
  out = train.fun(x,y) 
  fit = matrix(predict.fun(out,x),nrow=n)
  pred = matrix(predict.fun(out,x0),nrow=n0)
  m = ncol(pred)

  # Compute leave-one-out residuals with special function
  if (!is.null(special.fun) &&
      (is.null(mad.train.fun) && is.null(mad.predict.fun))) {
    res = abs(special.fun(x,y,out))
  }
  
  # Compute leave-one-out residuals by brute-force
  else {
    res = matrix(0,n,m)
    for (i in 1:n) {
      if (verbose) {
        cat(sprintf("\r%sProcessing training point %i (of %i) ...",txt,i,n))
        flush.console()
      }
     
      xx = x[-i,]; yy = y[-i]
      out.i = train.fun(xx,yy)
      res[i,] = abs(y[i] - predict.fun(out.i,x[i,,drop=F]))

      # Local scoring?
      if (!is.null(mad.train.fun) && !is.null(mad.predict.fun)) {
        r = abs(y[-i] - matrix(predict.fun(out.i,x[-i,,drop=F]),nrow=n-1))
        for (l in 1:m) {
          out.mad.i = mad.train.fun(x[-i,,drop=F],r[,l])
          res[i,l] = res[i,l] / mad.predict.fun(out.mad.i,x[i,,drop=F])
        }
      }
    }
    if (verbose) cat("\n")
  }

  # Compute quantiles and intervals
  lo = up = matrix(0,n0,m)

  for (l in 1:m) {
    # Local scoring?
    if (!is.null(mad.train.fun) && !is.null(mad.predict.fun)) {
      res.train = abs(y - fit[,l])
      mad.out = mad.train.fun(x,res.train)
      mad.x0 = mad.predict.fun(mad.out,x0)
    }
    else {
      mad.x0 = rep(1,n0)
    }
    
    q = quantile(res[,l],1-alpha)
    lo[,l] = pred[,l] - q * mad.x0
    up[,l] = pred[,l] + q * mad.x0
  }
  
  return(list(pred=pred,lo=lo,up=up,fit=fit))
}
