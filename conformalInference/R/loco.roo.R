#' Variable importance via excess test error
#'
#' Compute prediction intervals for the excess test error due to dropping a
#'   variable, measured in-sample and having valid in-sample coverage
#'
#' @param x Matrix of variables, of dimension (say) n x p.
#' @param y Vector of responses, of length (say) n.
#' @param train.fun A function to perform model training, i.e., to produce an
#'   estimator of E(Y|X), the conditional expectation of the response variable
#'   Y given features X. Its input arguments should be x: matrix of features,
#'   and y: vector of responses.
#' @param predict.fun A function to perform prediction for the (mean of the)
#'   responses at new feature values. Its input arguments should be out: output
#'   produced by train.fun, and newx: feature values at which we want to make
#'   predictions.
#' @param vars A list specifying the variables (indices between 1 and p) to
#'   for which variable importance should be investigated. Alternatively, if
#'   set equal to 0, the default, then all variables are investigated.
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
#'   define the first half of the data-split). Note that this split is common
#'   both to conformal.pred.roo (whether already run, or to be run internally,
#'   see below) and the excess error intervals computed by this function.
#'   Default is NULL, in which case the split is chosen randomly.
#' @param seed Integer to be passed to set.seed before defining the random 
#'   data-split to be used. Default is NULL, which effectively sets no seed.
#'   If both split and seed are passed, the former takes priority and the latter
#'   is ignored.
#' @param out.roo Output from running conformal.pred.roo on the given data set;
#'   default is NULL, which means that conformal.pred.roo will be run internally
#'   in order to compute in-sample prediction intervals. If out.roo is NULL,
#'   then the data-split common to conformal.pred.roo and this function will be
#'   determined by the split or seed arguments. If non-NULL, then the common
#'   data-split will be determined by that recorded in out.roo, and the split
#'   and seed arguments will be (silently) ignored.
#' @param verbose Should intermediate progress be printed out? Default is FALSE.
#'
#' @return A list with the following components: lo, up, vars, split,
#'   out.roo. The third component is a list of variables that were tested (i.e.,
#'   dropped one at a time, to compute excess error). The first two are arrays
#    of dimension n x length(vars) x m.
#'   The indices used for the half of the data-split are
#'   returned in split. Finally, for convenience, the output from running
#'   conformal.pred.roo is stored in out.roo.
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
#' @references "Distribution-Free Predictive Inference for Regression" by Lei,
#'   G'Sell, Rinaldo, Tibshirani, Wasserman (2018).
#' @example examples/ex.loco.roo.R
#' @export loco.roo

loco.roo = function(x, y, train.fun, predict.fun, vars=0,
  alpha=0.1, mad.train.fun=NULL, mad.predict.fun=NULL, split=NULL,
  seed=NULL, out.roo=NULL, verbose=FALSE) {

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

  inds = vector(mode="list",length=2)
  
  # Compute in-sample rank-one-out intervals ourselves
  if (is.null(out.roo)) {
    # If the user passed indices for the split, use them
    if (!is.null(split)) inds[[1]] = split
    # Otherwise make a random split
    else inds[[1]] = sample(1:n,floor(n/2))
  
    if (verbose) {
      cat(sprintf(paste("%sComputing in-sample rank-one-out conformal",
                        "prediction intervals ...\n"),txt))
    }

    out.roo = conformal.pred.roo(x,y,train.fun,predict.fun,alpha=alpha,
      mad.train.fun=mad.train.fun,mad.predict.fun=mad.predict.fun,
      split=inds[[1]],verbose=ifelse(verbose,paste0(txt,"\t"),FALSE))
  }
  # These have already been computed
  else {
    if (verbose) {
      cat(sprintf(paste("%sUsing given in-sample rank-one-out conformal",
                        "prediction intervals ...\n"),txt))
      cat(sprintf("%sSplitting data into parts of size %i and %i ...\n",txt,
                  floor(n/2),ceiling(n/2)))
    }
    inds[[1]] = out.roo$split
  }

  inds[[2]] = (1:n)[-inds[[1]]]
  pred = out.roo$pred
  m = ncol(pred)
  
  # Which variables to try dropping?
  if (length(vars) == 1 && vars == 0) vars = 1:p
  else {
    if (is.null(vars) || !is.numeric(vars) || min(vars) < 1 || 
        max(vars) > p || vars != round(vars)) {
      stop("vars must be a list of integers between 1 and p")
    }
    vars = unique(vars)
  }
  nv = length(vars)
  lo = up = array(0,dim=c(n,nv,m))
  
  # Train dropped variable models on one part, get residuals on the other
  for (k in 1:2) {
    i1 = inds[[k]]; n1 = length(i1)
    i2 = inds[[3-k]]; n2 = length(i2)

    cat(sprintf("%sHandling %s part of the data-split ...\n",
                txt,ifelse(k==1,"first","second")))
    
    for (j in 1:nv) {
      if (verbose) {
        cat(sprintf("\r%sDropping variable %i (%i of %i considered) ...",
                    txt,vars[j],j,nv))
        flush.console()
      }

      # Train on the first part, without variable 
      out.j = train.fun(x[i1,-vars[j],drop=F],y[i1])

      # Get predictions on the other, without variable
      pred.j = matrix(predict.fun(out.j,x[i2,-vars[j],drop=F]),nrow=n2)
      
      # Get the conformal intervals for differences in the residuals
      a = conformal.diff.int(pred[i2,,drop=F],pred.j,
        out.roo$lo[i2,,drop=F],out.roo$up[i2,,drop=F])
      lo[i2,j,] = a$lo
      up[i2,j,] = a$up 
    }
    
    if (verbose) cat("\n")
  }

  return(list(lo=lo,up=up,vars=vars,split=inds[[1]],out.roo=out.roo))
}

# Get the conformal interval for differences in residuals
# at mu2 and mu1, given a valid prediction interval [lo,up]
conformal.diff.int = function(mu1, mu2, lo, up) {
  n = nrow(mu1)
  m = ncol(mu1)
  v.lo = v.up = matrix(0,n,m)

  # Cases where mu1 is larger
  oo = mu2 <= mu1
  v.lo[oo] = trunc(lo,mu2,mu1)[oo]
  v.up[oo] = trunc(up,mu2,mu1)[oo]

  # Cases where mu2 is larger
  oo = mu1 <= mu2
  v.lo[oo] = -trunc(up,mu1,mu2)[oo]
  v.up[oo] = -trunc(lo,mu1,mu2)[oo]
  
  return(list(lo=v.lo,up=v.up))
}

###############################################################################

#' Variable importance (local) via excess test error, extended to right-
#' censored data
#' 
#' Compute prediction intervals for the excess test error due to dropping a
#'   variable, measured in-sample and having asymptotically valid in-sample 
#'   coverage. Extension to estimators of the restricted mean survival time with 
#'   right-censored data.
#'
#' @param x Matrix of variables, of dimension (say) n x p.
#' @param t Vector of responses (observed times), of length (say) n.
#' @param d Vector of responses (censoring indicator), of length (say) n.
#' @param tau Horizon of time.
#' @param train.fun A function to perform model training, i.e., to produce an
#'   estimator of E(min(T*,tau)|X), the conditional expectation of the true 
#'   event time T* given features X. Its input arguments should be x: matrix of 
#'   features, t: vector of observed times, d: vector of censoring indicators
#'   and tau: horizon of time. train.fun can be a function performing model 
#'   training for several (say m) regression models.
#' @param predict.fun A function to perform prediction for the (mean of the)
#'   responses at new feature values. Its input arguments should be out: output
#'   produced by train.fun, and newx: feature values at which we want to make
#'   predictions.
#' @param w Censoring weights. This should be a vector of length n. Default is 
#'   NULL, in which case censoring weights are computed using the censoring 
#'   model specified by the parameter cens.model, or by default with the 
#'   Kaplan-Meier model. 
#' @param cens.model Model used to estimate the censoring survival function.
#'   This should be equal to "km", "rsf" or "cox" (respectively Kaplan-Meier, 
#'   Random Survival Forests or Cox). Default is "km". 
#' @param vars A list specifying the variables (indices between 1 and p)
#'   for which variable importance should be investigated. Alternatively, if
#'   set equal to 0, the default, then all variables are investigated.
#' @param alpha Miscoverage level for the prediction intervals, i.e., intervals
#'   with coverage 1-alpha are formed. Default for alpha is 0.1.
#' @param mad.train.fun A function to perform training on the observed absolute 
#'   residuals R = |min(T,tau) - m(X)| and their associated censoring weights,  
#'   where m denotes the estimator produced by train.fun, to estimate E(R*|X)
#'   for the true absolute residual denoted R*.
#'   This is used to scale the conformal score, to produce a prediction interval
#'   with varying local width. The input arguments to mad.train.fun should be
#'   x: matrix of features, y: vector of observed absolute residuals, and w: 
#'   censoring weights. The default for mad.train.fun is NULL, which means that 
#'   no training is done on the absolute residuals, and the usual (unscaled) 
#'   conformal score is used. Note that if mad.train.fun is non-NULL, then so 
#'   must be mad.predict.fun (see next).
#' @param mad.predict.fun A function to perform prediction for the (mean of the)
#'   absolute residuals at new feature values. Its input arguments should be
#'   out: output produced by mad.train.fun, and newx: feature values at which we
#'   want to make predictions. The default for mad.predict.fun is NULL, which
#'   means that no local scaling is done for the conformal score, i.e., the
#'   usual (unscaled) conformal score is used.
#' @param split Indices that define the data-split to be used (i.e., the indices
#'   define the first half of the data-split). Note that this split is common
#'   both to conformal.pred.roo.surv (whether already run, or to be run inter-
#'   nally, see below) and the excess error intervals computed by this function.
#'   Default is NULL, in which case the split is chosen randomly.
#' @param seed Integer to be passed to set.seed before defining the random 
#'   data-split to be used. Default is NULL, which effectively sets no seed.
#'   If both split and seed are passed, the former takes priority and the latter
#'   is ignored.
#' @param out.roo.surv Output from running conformal.pred.roo.surv on the given  
#'   data set; default is NULL, which means that conformal.pred.roo.surv will be  
#'   run internally in order to compute in-sample prediction intervals. If  
#'   is NULL, then the data-split common to conformal.pred.roo.surv and this 
#'   out.roo.surv function will be determined by the split or seed arguments. If  
#'   non-NULL, then the common data-split will be determined by that recorded in  
#'   out.roo.surv, and the split and seed arguments will be (silently) ignored.
#' @param verbose Should intermediate progress be printed out? Default is FALSE.
#'
#' @return A list with the following components: lo, up, vars, split,
#'   out.roo. The third component is a list of variables that were tested (i.e.,
#'   dropped one at a time, to compute excess error). The first two are arrays
#'   of dimension n x length(vars) x m. The indices used for the half of the 
#'   data-split are returned in split. Finally, for convenience, the output 
#'   from running conformal.pred.roo.surv is stored in out.roo.
#'
#' @references "Distribution-Free Predictive Inference for Regression" by Lei,
#'   G'Sell, Rinaldo, Tibshirani, Wasserman (2018). See also "A Comprehensive 
#'   Framework for Evaluating Time to Event Predictions using the Restricted 
#'   Mean Survival Time" by Cwiling, Perduca, Bouaziz (2023) for the extension
#'   to right-censored data.
#' @export loco.roo.surv

loco.roo.surv = function(x, t, d, tau, train.fun, predict.fun, w=NULL, 
                         cens.model="km", vars=0, alpha=0.1, mad.train.fun=NULL, mad.predict.fun=NULL, 
                         split=NULL, seed=NULL, out.roo.surv=NULL, verbose=FALSE) {
  
  # Set up data
  x = as.matrix(x)
  t = as.numeric(t)
  d = as.numeric(d)
  if(!is.null(w)) w = as.numeric(w)
  n = nrow(x)
  p = ncol(x)
  
  # Check input arguments
  check.args.surv(x=x,t=t,d=d,tau=tau,x0=x,w=w,cens.model=cens.model,
                  alpha=alpha,train.fun=train.fun,predict.fun=predict.fun,
                  mad.train.fun=mad.train.fun,mad.predict.fun=mad.predict.fun)
  
  # Users may pass in a string for the verbose argument
  if (verbose == TRUE) txt = ""
  if (verbose != TRUE && verbose != FALSE) {
    txt = verbose
    verbose = TRUE
  }
  
  inds = vector(mode="list",length=2)
  
  # Compute in-sample rank-one-out intervals ourselves
  if (is.null(out.roo.surv)) {
    # If the user passed indices for the split, use them
    if (!is.null(split)) inds[[1]] = split
    # Otherwise make a random split
    else inds[[1]] = sample(1:n,floor(n/2))
    
    if (verbose) {
      cat(sprintf(paste("%sComputing in-sample rank-one-out conformal",
                        "prediction intervals ...\n"),txt))
    }
    
    out.roo.surv = conformal.pred.roo.surv(x,t,d,tau,train.fun,predict.fun,alpha=alpha,
                                           mad.train.fun=mad.train.fun,mad.predict.fun=mad.predict.fun, split=inds[[1]],
                                           w=w, cens.model=cens.model,verbose=ifelse(verbose,paste0(txt,"\t"),FALSE))
  }
  # These have already been computed
  else {
    if (verbose) {
      cat(sprintf(paste("%sUsing given in-sample rank-one-out conformal",
                        "prediction intervals ...\n"),txt))
      cat(sprintf("%sSplitting data into parts of size %i and %i ...\n",txt,
                  floor(n/2),ceiling(n/2)))
    }
    inds[[1]] = out.roo.surv$split
  }
  
  inds[[2]] = (1:n)[-inds[[1]]]
  pred = out.roo.surv$pred
  m = ncol(pred)
  
  # Which variables to try dropping?
  if (length(vars) == 1 && vars == 0) vars = 1:p
  else {
    if (is.null(vars) || !is.numeric(vars) || min(vars) < 1 || 
        max(vars) > p || vars != round(vars)) {
      stop("vars must be a list of integers between 1 and p")
    }
    vars = unique(vars)
  }
  nv = length(vars)
  lo = up = array(0,dim=c(n,nv,m))
  
  # Train dropped variable models on one part, get residuals on the other
  for (k in 1:2) {
    i1 = inds[[k]]; n1 = length(i1)
    i2 = inds[[3-k]]; n2 = length(i2)
    
    if (verbose) {
      cat(sprintf("%sHandling %s part of the data-split ...\n",
                  txt,ifelse(k==1,"first","second")))
    }
    
    for (j in 1:nv) {
      if (verbose) {
        cat(sprintf("\r%sDropping variable %i (%i of %i considered) ...",
                    txt,vars[j],j,nv))
        flush.console()
      }
      
      # Train on the first part, without variable 
      out.j = train.fun(x[i1,-vars[j],drop=F],t[i1],d[i1],tau)
      
      # Get predictions on the other, without variable
      pred.j = matrix(predict.fun(out.j,x[i2,-vars[j],drop=F],tau),nrow=n2)
      
      # Get the conformal intervals for differences in the residuals
      a = conformal.diff.int(pred[i2,,drop=F],pred.j,
                             out.roo.surv$lo[i2,,drop=F],out.roo.surv$up[i2,,drop=F])
      lo[i2,j,] = a$lo
      up[i2,j,] = a$up 
    }
    
    if (verbose) cat("\n")
  }
  
  return(list(lo=lo,up=up,vars=vars,split=inds[[1]],
              out.roo.surv=out.roo.surv))
}
