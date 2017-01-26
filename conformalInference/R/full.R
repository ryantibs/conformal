#' Conformal prediction intervals.
#'
#' Compute prediction intervals using conformal inference.
#'
#' @param x Matrix of features, of dimension (say) n x p.
#' @param y Vector of responses, of length (say) n.
#' @param x0 Matrix of features, each row being a point at which we want to
#'   form a prediction interval, of dimension (say) n0 x p. 
#' @param train.fun A function to perform model training, i.e., to produce an
#'   estimator of E(Y|X), the conditional expectation of the response variable
#'   Y given features X. Its input arguments should be x: matrix of features,
#'   y: vector of responses, and out: the output produced by a previous call
#'   to train.fun, at the \emph{same} features x. The function train.fun may
#'   (optionally) leverage this returned output for efficiency purposes. See
#'   details below.
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
#'   x: matrix of features, y: vector of absolute residuals, and out: the output
#'   produced by a previous call to mad.train.fun, at the \emph{same} features
#'   x. The function mad.train.fun may (optionally) leverage this returned
#'   output for efficiency purposes. See details below. The default for 
#'   mad.train.fun is NULL, which means that no training is done on the absolute
#'   residuals, and the usual (unscaled) conformal score is used. Note that if
#'   mad.train.fun is non-NULL, then so must be mad.predict.fun (next).
#' @param mad.predict.fun A function to perform prediction for the (mean of the)
#'   absolute residuals at new feature values. Its input arguments should be
#'   out: output produced by mad.train.fun, and newx: feature values at which we
#'   want to make predictions. The default for mad.predict.fun is NULL, which
#'   means that no local scaling is done for the conformal score, i.e., the
#'   usual (unscaled) conformal score is used.
#' @param num.grid.pts Number of grid points used when forming the conformal
#'   intervals (each grid point is a trial value for the interval). Default is
#'   100.
#' @param grid.factor Expansion factor used to define the grid for the conformal
#'   intervals. Default is 0.25, which means that the grid extends 25% beyond 
#'   the range of observed y's, on either side (to the left of the min, to the
#'   right of max). 
#' @param grid.method One of "linear", "floor", "ceiling". Default is "linear"
#'   representing the method used to interpolate the endpoints of the conformal
#'   conformal prediction interval from the grid of query values. The first
#'   "linear" corresponds to linear interpolation; the other options "floor"
#'   and "ceiling" choose the nearest and farthest grid endpoints, respectively,
#'   i.e., return the least conservative and most conservative prediction 
#'   intervals, respectively. When num.grid.pts is large, all three options will
#'   perform basically the same.
#' @param verbose Should intermediate progress be printed out? Default is FALSE.
#'
#' @return A list with the following components: pred, lo, up, fit. The first
#'   three are matrices of dimension m x n0, and the last is a matrix of
#'   dimension m x n. Recall that n0 is the number of rows of x0, and m is the
#'   number of tuning parameter values internal to predict.fun. In a sense, each
#'   of the m columns really corresponds to a different prediction function;
#'   see details below. Hence, the rows of the matrices pred, lo, up give
#'   the predicted value, and lower and upper confidence limits (from conformal
#'   inference), respectively, for the response at the n0 points given in
#'   x0. The rows of fit give the fitted values for the n points given in x.
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
#'   The third arugment to train.fun, as explained above, is the output produced
#'   by a previous call to train.fun, but importantly, at the \emph{same}
#'   features x. The function train.fun may (optionally) leverage this returned
#'   output for efficiency purposes. Here is how it is used, in this function: 
#'   when successively querying several trial values for the prediction
#'   interval, denoted, say, y0[j], j = 1,2,3,..., we set the out argument in
#'   train.fun to be the result of calling train.fun at the previous trial value
#'   y0[j-1]. Note of course that the function train.fun can also choose to 
#'   ignore the returned output out, and most default training functions made
#'   available in this package will do so, for simplicity. An exception is the
#'   training function produced by \code{\link{lm.funs}}. This will use this 
#'   output efficiently: the first time it trains by regressing onto a matrix of
#'   features x, it computes an appropriate Cholesky factorization; each
#'   successive time it is asked to train by regressing onto the same matrix x,
#'   it simply uses this Cholesky factorization (rather than recomputing one).
##   TODO implement and mention this for other functions too?
#'   Finally, the analogous explanation and discussion here applies to the out
#'   argument used by mad.train.fun.
#'   
#' @seealso \code{\link{conformal.pred.jack}},
#'   \code{\link{conformal.pred.split}}, \code{\link{conformal.pred.roo}}
#' @author Ryan Tibshirani
#' @references "Distribution-Free Predictive Inference for Regression" by 
#'   Jing Lei, Max G'Sell, Alessandro Rinaldo, Ryan Tibshirani, and Larry
#'   Wasserman, https://arxiv.org/pdf/1604.04173.pdf, 2016.
#' @example examples/ex.conformal.pred.R
#' @export conformal.pred

conformal.pred = function(x, y, x0, train.fun, predict.fun, alpha=0.1,
  mad.train.fun=NULL, mad.predict.fun=NULL, num.grid.pts=100, grid.factor=0.25,
  grid.method=c("linear","floor","ceiling"), verbose=FALSE) {

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
  if (length(num.grid.pts) != 1 || !is.numeric(num.grid.pts)
      || num.grid.pts <= 1 || num.grid.pts >= 1000
      || round(num.grid.pts) != num.grid.pts) {
    stop("num.grid.pts must be an integer between 1 and 1000")
  }
  check.pos.num(grid.factor)
  grid.method = match.arg(grid.method)
  
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
  
  # Trial values for y, empty lo, up matrices to fill
  rg = max(y)-min(y)
  yvals = seq(min(y)-grid.factor*rg, max(y)+grid.factor*rg,length=num.grid.pts)
  lo = up = matrix(0,n0,m)
  pvals = matrix(0,num.grid.pts,m)
  xx = rbind(x,rep(0,p))
  
  for (i in 1:n0) {
    if (verbose) {
      cat(sprintf("\r%sProcessing prediction point %i (of %i) ...",txt,i,n0))
      flush.console()
    }
    
    xx[n+1,] = x0[i,]
  
    for (j in 1:num.grid.pts) {
      yy = c(y,yvals[j])
      if (j==1) out = train.fun(xx,yy)
      else out = train.fun(xx,yy,out)
      r = abs(yy - matrix(predict.fun(out,xx),nrow=n+1))

      # Local scoring?
      if (!is.null(mad.train.fun) && !is.null(mad.predict.fun)) {
        for (l in 1:m) {
          if (j==1 && l==1) out.mad = mad.train.fun(xx,r[,l])
          else out.mad = mad.train.fun(xx,r[,l],out.mad)
          r[,l] = r[,l] / mad.predict.fun(out.mad,xx)
        }
      }

      rr = matrix(rep(r[n+1,],each=n+1),nrow=n+1)
      pvals[j,] = colMeans(rr <= r)
    }

    for (l in 1:m) {
      int = grid.interval(yvals,pvals[,l],alpha,right=T,method=grid.method)
      lo[i,l] = int$lo
      up[i,l] = int$up
    }
  }
  if (verbose) cat("\n")
  
  return(list(pred=pred,lo=lo,up=up,fit=fit))
}
