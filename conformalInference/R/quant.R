#' Full Conformal quantile regression prediction intervals.
#'
#' Compute prediction intervals using conformal quantile inference with Random Forest.
#'
#' @param x Matrix of features, of dimension (say) n x p.
#' @param y Vector of responses, of length (say) n.
#' @param x0 Matrix of features, each row being a point at which we want to
#'   form a prediction interval, of dimension (say) n0 x p. 
#' @param method Choose the method to compute the confidence intervals. 
#' The options are "classic" , without local scaling of the intervals, "median",
#' using the median instead of the pointwise evaluations (y), or "scaled", with
#' local scaling. Default is "scaled".
#' @param nthreads Set the number of threads to use (for parallel computation)
#' in the function quantregForest. Default is 8.
#' @param alpha Miscoverage level for the prediction intervals, i.e., intervals
#' with coverage 1-alpha are formed. Default for alpha is 0.1.
#' @param gamma Hyperparameter which defines the quantile levels used in the
#' method. Default is alpha/2.
#' @param split Indices that define the data-split to be used (i.e., the indices
#'   define the first half of the data-split, on which the model is trained).
#'   Default is NULL, in which case the split is chosen randomly.
#' @param seed Integer to be passed to set.seed before defining the random 
#'   data-split to be used. Default is NULL, which effectively sets no seed.
#'   If both split and seed are passed, the former takes priority and the latter
#'   is ignored.
#' @param verbose Should intermediate progress be printed out? Default is FALSE.
#' @param w Weights, in the case of covariate shift. This should be a vector of
#'   length n+n0, giving the weights (i.e., ratio of test to training feature
#'   densities), at each of n+n0 the training and test points. Default is NULL,
#'   which means that we take all weights to be 1.
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
#'   intervals, i.e., the grid points are taken to be equally spaced in between
#'   -grid.factor*max(abs(y)) and grid.factor*max(abs(y)). Default is 1.25. In
#'   this case (and with exchangeable data, thus unity weights) the restriction
#'   of the trial values to this range costs at most 1/(n+1) in coverage. See
#'   details below.
#'
#' @return A list with the following components: lo, up,i.e.two 
#'  vectors of length n0, while split contains the training indices.
#' 
#' @details This function is based on the package \code{\link{quantregForest}} to
#'  compute the Quantile Regression. 
#'  If this package is not installed, then the function will abort.
#'  
#' @details Since the computation can be quite long, we exploited the package
#' \code{\link{pbapply}} to provide progress bar for loops. If this package 
#' is not installed, then the classical apply function will be used.
#'
#' @details The function may create invalid confidence intervals 
#' (this problem is described in the Sesia, Candès reference paper and occurs when the
#' selected method is "classic"). If this happens, an error is displayed suggesting 
#' to fine-tune the value of the hyperparameter gamma or to change the input method.
#'
#' @references The quantile nonconformity measures are taken fom  "A comparison of 
#' some conformal quantile regression methods" 
#' by Sesia, Candés (2020).
#' Moreover the "classic","median" and "scaled" methods are respectively taken from
#' "Conformalized quantile regression" by Romano, Patterson, Candès (2019), from the
#' previous Sesia, Candès paper and from "A (tight) upper bound for the length of
#' confidence intervals with conditional coverage" by Kivaranovic, Leeb (2020).
#'   
#' @export conformal.quant


conformal.quant = function(x, y, x0, method="scaled", nthreads = 8,
                           alpha=0.1, 
                           gamma=alpha/2, verbose=FALSE,w=NULL,
                           mad.train.fun=NULL,
                           mad.predict.fun=NULL, num.grid.pts=25, grid.factor=1.25) {
  
  
  
  # Set up data
  x = as.matrix(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)
  x0 = matrix(x0,ncol=p)
  n0 = nrow(x0)
  eps = 1e-6
  what=c(gamma,1-gamma)
  
  check.quant(x=x, y=y, x0=x0, method=method, nthreads, alpha=alpha, gamma=gamma, w=w,
              mad.train.fun=mad.train.fun,mad.predict.fun+mad.predict.fun,
              num.grid.pts=num.grid.pts, grid.factor=grid.factor)  
  
  # Check the weights
  if (is.null(w)) w = rep(1,n+n0)
  
  # Users may pass in a string for the verbose argument
  if (verbose == TRUE) txt = ""
  if (verbose != TRUE && verbose != FALSE) {
    txt = verbose
    verbose = TRUE
  }
  
  if (verbose) cat(sprintf("%sInitial training on full data set ...\n",txt))
  
  
  #create grid of values for y
  ymax = max(abs(y))
  yvals = seq(-grid.factor*ymax, grid.factor*ymax,length=num.grid.pts)
  lo = up = rep(0,n0)
  #qvals = rvals = rep(0,num.grid.pts)
  xx = rbind(x,rep(0,p))
  
  #compute scores and p-values
  outer<-t(pbar.sapply(1:n0, function(i) {
    if (verbose) {
      cat(sprintf("\r%sProcessing prediction point %i (of %i) ...",txt,i,n0))
      flush.console()
    }
    
    xx[n+1,] = x0[i,]
    ww = c(w[1:n],w[n+i])
    
    # Refit for each point in yvals, compute conformal p-value
    inner<-t(sapply (1:num.grid.pts, function(j) {
      yy = c(y,yvals[j])
      out=quantregForest::quantregForest(xx,yy, nthreads)
      
      #Choose quantile regression method
      switch (method,
              
              "classic" ={ 
                fit = matrix(predict(out,xx,what), nrow=n+1)
                band = fit %*% diag(c(1,-1)) + cbind(-yy, yy)
                r=apply(band,1,max)
              },
              
              "median" = {
                what_med = c(what[1],0.5,what[2])
                fit = matrix(predict(out,xx,what_med), nrow=n+1)
                error_low = (fit[,1] - yy)/(fit[,2] - fit[,1] +eps )
                error_high = (yy - fit[,3])/(fit[,1] - fit[,2] +eps )
                band = cbind(error_low,error_high)
                r=apply(band,1,max)
              },
              
              "scaled" = {
                fit = matrix(predict(out,xx,what), nrow=n+1)
                scaling_factor = fit[,2] - fit[,1] + eps
                error_low = (fit[,1] - yy)/scaling_factor
                error_high = (yy - fit[,2])/scaling_factor
                band = cbind(error_low,error_high)
                r=apply(band,1,max)
              },
      )
      
      
      # Local scoring?
      if (!is.null(mad.train.fun) && !is.null(mad.predict.fun)) {
        if (j==1 ) out.mad = mad.train.fun(xx,r)
        else out.mad = mad.train.fun(xx,r,out.mad)
        r = r / mad.predict.fun(out.mad,xx)
        
      }
      
      return(c(weighted.quantile(r,prob=1-alpha,w=ww),r[n+1]))
    }))
    
    qvals=inner[,1]
    rvals=inner[,2]
    
    int = grid.interval(yvals,rvals,qvals)
    return(c(int$lo,int$up))
    
  }))
  
  lo<-outer[,1]
  up<-outer[,2]
  
  if (verbose) cat("\n")
  
  
  if(sum(up - lo <0)>0)
    stop("Some confidence intervals are not valid. 
    Try to tune the hyperparameter gamma 
         or to change the value of the input parameter method\n")
  
  return(list(lo=lo,up=up))
}