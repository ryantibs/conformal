#' Split Conformal quantile regression prediction intervals.
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
#'
#' @return A list with the following components: lo, up, split. The first two 
#' are vectors of length n0, while split contains the training indices.
#' 
#' @details This function is based on the package \code{\link{quantregForest}} to
#'  compute the Quantile Regression. 
#'  If this package is not installed, then the function will abort. 
#'
#' @details The function may create invalid confidence intervals 
#' (this problem is described in the Sesia, Candès reference paper and occurs when the
#' selected method is "classic"). If this happens, an error is displayed suggesting 
#' to fine-tune the value of the hyperparameter gamma or to change the input method.
#'
#' @references See "A comparison of some conformal quantile regression methods" 
#' by Sesia, Candés (2020), which also links to the code implementation in python.
#' Moreover the "classic","median" and "scaled" methods are respectively taken from
#' "Conformalized quantile regression" by Romano, Patterson, Candès (2019), from the
#' previous Sesia, Candès paper and from "A (tight) upper bound for the length of
#' confidence intervals with conditional coverage" by Kivaranovic, Leeb (2020).
#'   
#' @export conformal.quant.split


conformal.quant.split = function(x, y, x0, method="scaled", nthreads = 8,
                                 alpha=0.1, 
                                 gamma=alpha/2, split=NULL, 
                                 seed=NULL, verbose=FALSE) {
  
  
  possible_functions=c("classic","median","scaled")
  if (is.null(method) || method %in% possible_functions==FALSE) {
    stop(c("The 'method' argument is not correct. 
             Please select one of the following:",
           paste(possible_functions,collapse=", "),"."))
  }
  
  
  # Set up data
  x = as.matrix(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)
  x0 = matrix(x0,ncol=p)
  n0 = nrow(x0)
  eps = 1e-6
  what=c(gamma,1-gamma)
  
  
  # Check input arguments
  if (is.null(x) || !is.numeric(x)) stop("x must be a numeric matrix")
  if (nrow(x) != length(y)) stop("nrow(x) and length(y) must match")
  if (is.null(x0) || !is.numeric(x0)) stop("x0 must be a numeric matrix")
  check.num.01(alpha)
  check.num.01(gamma)
  
  
  
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
  
  out = quantregForest::quantregForest(x[i1,],y[i1], nthreads)
  
  
  # Train quantiles on first part
  if (verbose) {
    cat(sprintf("%sComputing residuals and quantiles on second part ...\n",txt))
  }
  
  
  # Switch among the various methods
  switch (method,
          
          "classic" ={ 
            pred = matrix(predict(out,x0,what), nrow=nrow(x0))
            fit = matrix(predict(out,x[i2,],what), nrow=length(i2))
            
            band = fit %*% diag(c(1,-1)) + cbind(-y[i2], y[i2])
            r=apply(band,1,max)
            q = quantile(c(r,Inf),probs = 1-alpha)
            lo=pred[,1]- q 
            up=pred[,2]+ q
          },
          
          "median" = {
            what_med = c(what[1],0.5,what[2])
            pred = matrix(predict(out,x0,what_med), nrow=nrow(x0))
            fit = matrix(predict(out,x[i2,],what_med), nrow=length(i2))
            
            error_low = (fit[,1] - y[i2])/(fit[,2] - fit[,1] +eps )
            error_high = (y[i2] - fit[,3])/(fit[,1] - fit[,2] +eps )
            band = cbind(error_low,error_high)
            r=apply(band,1,max)
            q = quantile(c(r,Inf),probs = 1-alpha)
            lo=pred[,1]- q * (pred[,2] - pred[,1] + eps)
            up=pred[,2]+ q * (pred[,3] - pred[,2] + eps)
          },
          
          "scaled" = {
            pred = matrix(predict(out,x0,what), nrow=nrow(x0))
            fit = matrix(predict(out,x[i2,],what), nrow=length(i2))
            
            scaling_factor = fit[,2] - fit[,1] + eps
            error_low = (fit[,1] - y[i2])/scaling_factor
            error_high = (y[i2] - fit[,2])/scaling_factor
            band = cbind(error_low,error_high)
            r=apply(band,1,max)
            q = quantile(c(r,Inf),probs = 1-alpha)
            pred_factor = pred[,2] - pred[,1] + eps
            lo=pred[,1]- q * pred_factor
            up=pred[,2]+ q * pred_factor
          },
  )
  
  
  
  if(sum(up - lo <0)>0)
    stop("Some confidence intervals are not valid. 
    Try to tune the hyperparameter gamma 
         or to change the value of the input parameter method\n")
  
  return(list(lo=lo,up=up,split=i1))
}