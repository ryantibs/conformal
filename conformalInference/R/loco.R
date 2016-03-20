#' Feature importance via sample splitting.
#'
#' Compute confidence intervals for mean (or median) excess prediction error
#'   after leaving out one feature.
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
#' @param active.fun A function to extract the active features for a 
#'   
#' @param alpha Miscoverage level for the confidence intervals, i.e., intervals
#'   with coverage 1-alpha are formed. Default for alpha is 0.1.
#' @param split Indices that define the data-split to be used (i.e., the indices
#'   define the first half of the data-split, on which the model is trained).
#'   Default is NULL, in which case the split is chosen randomly.
#' @param seed Integer to be passed to set.seed before defining the random 
#'   data-split to be used. Default is NULL. If both split and seed are passed,
#'    the former takes priority and the latter is ignored.
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
#' @author Ryan Tibshirani, Larry Wasserman, and CMU Conformal Inference Team
#' @references \url{http://www.stat.cmu.edu}
#' @example examples/ex.loco.R
#' @export loco

loco = function(x, y, train.fun, predict.fun, active.fun, alpha=0.1,
  split=NULL, seed=NULL, verbose=FALSE) {

  # Set up data
  x = as.matrix(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)
  
  # Check input arguments
  check.args(x=x,y=y,x0=x,alpha=alpha,train.fun=train.fun,
             predict.fun=predict.fun)

  # Users may pass in a string for the verbose argument
  if (verbose == TRUE) txt = ""
  if (verbose != TRUE && verbose != FALSE) {
    txt = verbose
    verbose = TRUE
  }
  
  # Users may pass in a string for the verbose argument
  if (verbose == TRUE) txt = ""
  if (verbose != TRUE && verbose != FALSE) {
    txt = verbose
    verbose = TRUE
  }
  
  # If the user passed indices for the split, use them
  if (!is.null(split)) i1 = split
  # Otherwise make a random split
  else { set.seed(seed); i1 = sample(1:n,floor(n/2)) }
  i2 = (1:n)[-i1]
  n1 = length(i1)
  n2 = length(i2)
  
  if (verbose) {
    cat(sprintf("%sSplitting data into parts of size %i and %i ...\n",
                txt,n1,n2))
    cat(sprintf("%sInitial training on first part ...\n",txt))
  }
  
  # Train on first part
  out = train.fun(x[i1,,drop=F],y[i1])
  active = active.fun(out)
  master = unique(unlist(active))
  m = length(active)
  J = length(master)
  
  if (verbose) {
    cat(sprintf("%Initial residuals on second part ...\n",txt))
  }
  
  # Get residuals on second
  res = abs(y[i2] - matrix(predict.fun(out,x[i2,,drop=F]),nrow=n2))
  res.drop = vector(mode="list",length=J)

  # Re-fit after dropping each feature in the master list
  for (j in 1:J) {
    if (verbose) {
      cat(sprintf(paste("%s\tRe-fitting after dropping feature %i (%i of %i",
                        "considered) ...\n"),txt,master[j],j,J))
      flush.console()
    }

    # Train on the first part, without variable 
    out.j = train.fun(x[i1,-active[j],drop=F],y[i1])

    # Get predictions on the other, without variable
    res.drop[[j]] = abs(y[i2] - matrix(predict.fun(out.j,x[i2,-active[j],
              drop=F]),nrow=n2))  
  }
  if (verbose) cat("\n")

  # Compute p-values and confidence intervals at each tuning step
  inf.mean = inf.sign = inf.wilcox = vector(mode="list",length=m)
  
  for (l in 1:m) {
    if (verbose) {
      cat(sprintf("%s\tPerforming LOCO inference at tuning step %i (of %i)",
                  txt,l,m))
      flush.console()
    }

    k = length(active[[l]])
    inf.mean[[l]] = inf.sign[[l]] = inf.wilcox[[l]] = matrix(0,k,3)
    for (j in 1:k) {
      j.master = which(master == active[j])
      z = res.drop[[j.master]][,l] - res[,l]
      inf.mean[[l]][j,] = my.z.test(z,alpha)
      inf.sign[[l]][j,] = my.sign.test(z,alpha)
      inf.wilcox[[l]][j,] = my.wilcox.test(z,alpha)
    }
  }
  if (verbose) cat("\n")

  return(list(inf.mean=inf.mean,inf.sign=inf.sign,inf.wilcox=inf.wilcox))
}

# Mean inference: One-sided p-value but a two-sided confidence interval
my.z.test = function(z, alpha){
  n = length(z)
  s = sd(z)
  m = mean(z)
  pval = 1-pnorm(m/s*sqrt(n))
  q = qnorm(1-alpha/2)
  left  = m - q*s/sqrt(n)
  right = m + q*s/sqrt(n)
  return(c(pval,left,right))
}

# Median inference: one-sided p-value but a two-sided confidence interval
my.sign.test = function(z, alpha){
  n = length(z)
  s = sum(z>0)
  pval = 1-pbinom(s-1,n,0.5)
  k = n - qbinom(1-alpha/2,n,0.5) # Want P(B <= k) <= alpha/2
  zz  = sort(z)
  left  = zz[k]
  right = zz[n-k+1]
  return(c(pval,left,right))
}

# Median inference: one-sided p-value but a two-sided confidence interval
my.wilcox.test = function(z, alpha){
  pval = wilcox.test(z, alternative="greater")$p.value
  out = wilcox.test(z, conf.int=TRUE, conf.level=1-alpha)
  left = out$conf.int[1]
  right = out$conf.int[2]
  return(c(pval,left,right))
}




