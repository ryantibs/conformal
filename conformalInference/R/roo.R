#' In-sample split conformal prediction intervals. 
#'
#' Compute prediction intervals, having valid in-sample coverage, using
#'   rank-one-out (ROO) split conformal inference. 
#'
#' @param x Matrix of features, of dimension (say) n x p.
#' @param y Vector of responses, of length (say) n.
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
#'   define the first half of the data-split). Default is NULL, in which case
#'   the split is chosen randomly.
#' @param seed Integer to be passed to set.seed before defining the random 
#'   data-split to be used. Default is NULL, which effectively sets no seed.
#'   If both split and seed are passed, the former takes priority and the latter
#'   is ignored.
#' @param verbose Should intermediate progress be printed out? Default is FALSE.
#'
#' @return A list with the following components: pred, lo, up, fit, split,
#'   fit.all, out.all.
#'   The first four are matrices of dimension n x m. In a sense, each
#'   of the m columns really corresponds to a different prediction function; 
#'   see details below. Hence, the rows of the matrices pred, lo, up give the
#'   predicted value, and lower and upper confidence limits (from rank-one-out
#'   split conformal inference), respectively, for the response at the n points
#'   given in x. The rows of the matrix fit give the fitted values for the n
#'   points given in x. The indices used for the half of the data-split are
#'   returned in split. Finally, for convenience, the output from running
#'   train.out on the entire (unsplit) data set x,y is stored in out.all, and
#'   fit.all is an n x m matrix whose rows store the fitted values from out.all.
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
#'   The rank-one-out or ROO variant of split conformal prediction produces
#'   intervals with a special in-sample average coverage property. Namely, if
#'   we were to observe new
#'   response values y* at the given feature values x, then the ROO method
#'   produces a prediction band \eqn{\hat{C}(x)} with the asymptotic in-sample
#'   average coverage property
#'   \deqn{\liminf_{n \to \infty} \frac{1}{n} \sum_{i=1}^n 
#'   P \big(Y_i^* \in \hat{C}(X_i) \big) \geq 1-\alpha.}
#'   The usual split conformal band would not necessarily share this in-sample
#'   property. Of course, the split conformal band has valid predictive (i.e.,
#'   out-of-sample) coverage, and the ROO split conformal variant maintains
#'   this property as well.
#'
#' @seealso \code{\link{conformal.pred}}, 
#'   \code{\link{conformal.pred.jack}}, \code{\link{conformal.pred.split}}
#' @author Ryan Tibshirani
#' @references "Distribution-Free Predictive Inference for Regression" by 
#'   Jing Lei, Max G'Sell, Alessandro Rinaldo, Ryan Tibshirani, and Larry
#'   Wasserman, https://arxiv.org/pdf/1604.04173.pdf, 2016.
#' @example examples/ex.conformal.pred.roo.R
#' @export conformal.pred.roo

conformal.pred.roo = function(x, y, train.fun, predict.fun, alpha=0.1,
  mad.train.fun=NULL, mad.predict.fun=NULL, split=NULL, seed=NULL,
  verbose=FALSE) {

  # Set up data
  x = as.matrix(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)
  
  # Check input arguments
  check.args(x=x,y=y,x0=x,alpha=alpha,train.fun=train.fun,
             predict.fun=predict.fun,mad.train.fun=mad.train.fun,
             mad.predict.fun=mad.predict.fun)

  # Users may pass in a string for the verbose argument
  if (verbose == TRUE) txt = ""
  if (verbose != TRUE && verbose != FALSE) {
    txt = verbose
    verbose = TRUE
  }
  
  if (verbose) cat(sprintf("%sInitial training on full data set ...\n",txt))
  
  # Train and fit on full data set
  out.all = train.fun(x,y) 
  fit.all = matrix(predict.fun(out.all,x),nrow=n)
  m = ncol(fit.all)
  pred = lo = up = fit = matrix(0,n,m)
  
  # If the user passed indices for the split, use them
  inds = vector(mode="list",length=2)
  if (!is.null(split)) inds[[1]] = split
  # Otherwise make a random split
  else {
    if (!is.null(seed)) set.seed(seed)
    inds[[1]] = sample(1:n,floor(n/2))
  }
  inds[[2]] = (1:n)[-inds[[1]]]
  
  if (verbose) {
    cat(sprintf("%sSplitting data into parts of size %i and %i ...\n",txt,
                floor(n/2),ceiling(n/2)))
  }

  # Train one split, and get residuals on the other
  for (k in 1:2) {
    i1 = inds[[k]]; n1 = length(i1)
    i2 = inds[[3-k]]; n2 = length(i2)

    if (verbose) {
      cat(sprintf("%sTraining on %s part ...\n",txt,
                  ifelse(k==1,"first","second")))
    }

    # Train on first part
    out = train.fun(x[i1,,drop=F],y[i1])
    yhat1 = matrix(predict.fun(out,x[i1,,drop=F]),nrow=n1)
    yhat2 = matrix(predict.fun(out,x[i2,,drop=F]),nrow=n2)
    fit[i1,] = yhat1
    pred[i2,] = yhat2

    if (verbose) {
      cat(sprintf("%sComputing residuals and quantiles on %s part ...\n",
                  txt,ifelse(k==1,"first","second")))
    }
    
    # Get residuals and quantiles on second
    res = abs(y[i2] - yhat2)
    qmat = matrix(0,n2,m)
    
    for (l in 1:m) {
      # Local scoring?
      if (!is.null(mad.train.fun) && !is.null(mad.predict.fun)) {
        res.train = abs(y[i1] - yhat1[,l])
        mad.out = mad.train.fun(x[i1,,drop=F],res.train)
        mad.x2 = mad.predict.fun(mad.out,x[i2,,drop=F])
        res = res[,l] / mad.x2
      }
      else {
        mad.x2 = rep(1,n2)
      }
      
      # Sort once, update cleverly
      o = order(res[,l])
      sorted.res = res[o,l]
      oo = order(o)
      
      for (i in 1:n2) {
        qmat[i,l] = conformal.quantile(sorted.res[-oo[i]],alpha) * mad.x2[i]
      }
    }
    
    lo[i2,] = yhat2 - qmat
    up[i2,] = yhat2 + qmat
  }
  
  return(list(pred=pred,lo=lo,up=up,fit=fit,split=inds[[1]],
              fit.all=fit.all,out.all=out.all))
}
