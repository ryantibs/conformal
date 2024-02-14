#' Variable importance via mean (or median) excess test error
#'
#' Compute confidence intervals for mean (or median) excess test error due to 
#'   dropping a variable
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
#' @param active.fun A function which takes the output of train.fun, and reports
#'   which features are active for each fitted model contained in this output.
#'   Its only input argument should be out: output produced by train.fun.
#' @param alpha Miscoverage level for the confidence intervals, i.e., intervals
#'   with coverage 1-alpha are formed. Default for alpha is 0.1.
#' @param bonf.correct Should a Bonferroni correction be applied to the p-values
#'   and confidence intervals? Default is TRUE.
#' @param split Indices that define the data-split to be used (i.e., the indices
#'   define the first half of the data-split, on which the model is trained).
#'   Default is NULL, in which case the split is chosen randomly.
#' @param seed Integer to be passed to set.seed before defining the random 
#'   data-split to be used. Default is NULL, which effectively sets no seed.
#'   If both split and seed are passed, the former takes priority and the latter
#'   is ignored.
#' @param verbose Should intermediate progress be printed out? Default is FALSE.
#'
#' @return A list with the following components: inf.z, inf.sign, inf.wilcox,
#'   active, master, bonf.correct. The first three are lists, containing the
#'   results of LOCO inference with the Z-test, sign text, and Wilcoxon signed
#'   rank test, respectively. More details on these tests are given below. These
#'   lists have one element per tuning step inherent to the training and
#'   prediction functions, train.fun and predict.fun. The fourth returned
#'   component active is a list, with one element per tuning step, that reports
#'   which features are active in the corresponding fitted model. The fifth
#'   returned component master collects all active features across all tuning
#'   steps, for easy reference. The last component signals whether a Bonferroni
#'   correction has been applied.
#'
#' @details In leave-one-covariate-out or LOCO inference, the training data is
#'   split in two parts, and the first part is used to train a model, or some
#'   number of models across multiple tuning steps (e.g., indexed by different 
#'   tuning parameter values lambda in the lasso, or different steps along the
#'   forward stepwise path). For each model, each variable is left out one at
#'   time from the first part of the data, and the entire training procedure is
#'   repeated without this variable in consideration. Residuals are computed on
#'   both from the original fitted model, and the refitted model without the
#'   variable in consideration. A pairwise difference between the latter and
#'   former residuals is computed, and either a Z-test, sign test, or Wilcoxon
#'   signed rank test is performed to test either the mean or median difference
#'   here being zero. The Z-test is of course approximately valid under the
#'   conditions needed for the CLT; the sign test is distribution-free and exact
#'   in finite samples (only assumes continuity of the underlying distribution);
#'   the Wilcoxon signed rank test is also distribution-free and exact in finite
#'   samples, but may offer more power (it assumes both continuity and symmetry
#'   of the underlying distribution).
#'
#'   A couple other important notes: p-values here are from a one-sided test of
#'   the target parameter (mean or median excess test error) being equal to zero
#'   versus greater than zero. Confidence intervals are from inverting the
#'   two-sided version of this test.
#'
#' @references "Distribution-Free Predictive Inference for Regression" by Lei,
#'   G'Sell, Rinaldo, Tibshirani, Wasserman (2018).
#' @example examples/ex.loco.R
#' @export loco

loco = function(x, y, train.fun, predict.fun, active.fun, alpha=0.1,
                bonf.correct=TRUE, split=NULL, seed=NULL, verbose=FALSE) {
  
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
    cat(sprintf("%sInitial training on first part ...\n",txt))
  }
  
  # Train on first part
  out = train.fun(x[i1,,drop=F],y[i1])
  active = active.fun(out)
  master = unique(unlist(active))
  m = length(active)
  J = length(master)
  
  if (verbose) {
    cat(sprintf("%sInitial residuals on second part ...\n",txt))
  }
  
  # Get residuals on second
  res = abs(y[i2] - matrix(predict.fun(out,x[i2,,drop=F]),nrow=n2))
  res.drop = vector(mode="list",length=J)
  
  # Re-fit after dropping each feature in the master list
  for (j in Seq(1,J)) {
    if (verbose) {
      cat(sprintf(paste("\r%sRe-fitting after dropping feature %i (%i of %i",
                        "considered) ..."),txt,master[j],j,J))
      flush.console()
    }
    
    # Train on the first part, without variable 
    out.j = train.fun(x[i1,-master[j],drop=F],y[i1])
    
    # Get predictions on the other, without variable
    res.drop[[j]] = abs(y[i2] - matrix(predict.fun(out.j,x[i2,-master[j],
                                                           drop=F]),nrow=n2))  
  }
  if (verbose) cat("\n")
  
  # Compute p-values and confidence intervals at each tuning step
  inf.z = inf.sign = inf.wilcox = vector(mode="list",length=m)
  
  for (l in 1:m) {
    if (verbose) {
      if (m>1) {
        cat(sprintf(paste("\r%sPerforming LOCO analysis at tuning step %i",
                          "(of %i) ..."),txt,l,m))
        flush.console()
      }
      else cat(sprintf("%sPerforming LOCO analysis ...",txt))
    }
    
    k = length(active[[l]])
    inf.z[[l]] = inf.sign[[l]] = inf.wilcox[[l]] = matrix(0,k,3)
    for (j in Seq(1,k)) {
      j.master = which(master == active[[l]][j])
      z = res.drop[[j.master]][,l] - res[,l]
      inf.z[[l]][j,] = my.z.test(z,alpha,k,bonf.correct)
      inf.sign[[l]][j,] = my.sign.test(z,alpha,k,bonf.correct)
      inf.wilcox[[l]][j,] = my.wilcox.test(z,alpha,k,bonf.correct)
    }
    
    rownames(inf.z[[l]]) = rownames(inf.sign[[l]]) =
      rownames(inf.wilcox[[l]]) = active[[l]]
    colnames(inf.z[[l]]) = colnames(inf.sign[[l]]) =
      colnames(inf.wilcox[[l]]) = c("P-value", "LowConfPt", "UpConfPt")
  }
  if (verbose) cat("\n")
  
  out = list(inf.z=inf.z,inf.sign=inf.sign,inf.wilcox=inf.wilcox,
             active=active, master=master, bonf.correct=bonf.correct)
  class(out) = "loco"
  return(out)
}

#' Print function for loco object.
#'
#' Summary and print the results of a set of simulations, stored an object
#' of class sim (produced by sim.master).
#'
#' @param x The loco object.
#' @param test One of "z", "sign", "wilcox", or "all", describing which test 
#'   results to show. If "z", then the results from only Z-tests are displayed;
#'   if "sign", then only sign tests; if "wilcox", then only Wilcoxon tests. If
#'   "all", then the results from all tests are displayed. Default is "wilcox".
#' @param digits Number of digits to display. Default is 3. 
#' @param ... Other arguments (currently not used).
#'
#' @export 

print.loco = function(x, test=c("wilcox","sign","z","all"), digits=3, ...) {
  test = match.arg(test)
  
  cat(paste("\nDisplaying results of LOCO analysis. Notes:\n",
            "\n- P-values are from a one-sided test of the target parameter",
            "\n  (mean or median excess test error) being equal to zero versus",
            "\n  greater than zero.",
            "\n- Confidence intervals are from inverting the two-sided version",
            "\n  of this test.",
            "\n- Bonferonni correction:",
            ifelse(x$bonf.correct, "yes.", "no."), "\n"))
  
  m = length(x$inf.z)
  for (i in 1:m) {
    if (m>1) cat(sprintf("\n------- Tuning step %i -------\n",i))
    
    if (test=="z" || test=="all") {
      cat("\nZ-test results:\n\n")
      tab = round(x$inf.z[[i]],digits=digits)
      code = get.signif.code(x$inf.z[[i]][,1])
      tab = cbind(rownames(tab),tab,code)
      colnames(tab)[1] = "Var"
      colnames(tab)[5] = ""
      rownames(tab) = rep("",nrow(tab))
      print(tab,quote=FALSE)
    }
    
    if (test=="sign" || test=="all") {
      cat("\nSign test results:\n\n")
      tab = round(x$inf.sign[[i]],digits=digits)
      code = get.signif.code(x$inf.sign[[i]][,1])
      tab = cbind(rownames(tab),tab,code)
      colnames(tab)[1] = "Var"
      colnames(tab)[5] = ""
      rownames(tab) = rep("",nrow(tab))
      print(tab,quote=FALSE)
    }
    
    if (test=="wilcox" || test=="all") {
      cat("\nWilcoxon test results:\n\n")
      tab = round(x$inf.wilcox[[i]],digits=digits)
      code = get.signif.code(x$inf.wilcox[[i]][,1])
      tab = cbind(rownames(tab),tab,code)
      colnames(tab)[1] = "Var"
      colnames(tab)[5] = ""
      rownames(tab) = rep("",nrow(tab))
      print(tab,quote=FALSE)
    }
    
    cat("\nSignificance codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  }
}

# Mean inference: One-sided p-value but a two-sided confidence interval
my.z.test = function(z, alpha, k, bonf.correct=TRUE){
  n = length(z)
  s = sd(z)
  m = mean(z)
  pval = 1-pnorm(m/s*sqrt(n))
  
  # Apply Bonferroni correction for k tests
  if (bonf.correct) {
    pval = min(k*pval,1)
    alpha = alpha/k
  }
  
  q = qnorm(1-alpha/2)
  left  = m - q*s/sqrt(n)
  right = m + q*s/sqrt(n)
  return(c(pval,left,right))
}

# Median inference: one-sided p-value but a two-sided confidence interval
my.sign.test = function(z, alpha, k, bonf.correct=TRUE){
  n = length(z)
  s = sum(z>0)
  pval = 1-pbinom(s-1,n,0.5)
  
  # Apply Bonferroni correction for k tests
  if (bonf.correct) {
    pval = min(k*pval,1)
    alpha = alpha/k
  }
  
  j = n - qbinom(1-alpha/2,n,0.5) # Want P(B <= j) <= alpha/2
  zz  = sort(z)
  left  = zz[j]
  right = zz[n-j+1]
  return(c(pval,left,right))
}

# Median inference: one-sided p-value but a two-sided confidence interval
my.wilcox.test = function(z, alpha, k, bonf.correct=TRUE){
  pval = wilcox.test(z, alternative="greater", exact=TRUE)$p.value
  
  # Apply Bonferroni correction for k tests
  if (bonf.correct) {
    pval = min(k*pval,1)
    alpha = alpha/k
  }
  
  out = wilcox.test(z, conf.int=TRUE, conf.level=1-alpha, exact=TRUE)
  left = out$conf.int[1]
  right = out$conf.int[2]
  return(c(pval,left,right))
}

# Significance code function
get.signif.code = function(v) {
  code = rep("",length(v))
  code[v < 0.1] = "."
  code[v < 0.05] = "*"
  code[v < 0.01] = "**"
  code[v < 0.001] = "***"
  return(code)
}

###############################################################################

#' Variable importance (global) via median excess test error extended to right-
#' censored data
#'
#' Compute confidence intervals for median excess test error due to 
#'   dropping a variable. Extension to estimators of the restricted mean
#'   survival time with right-censored data.
#'
#' @param x Matrix of features, of dimension (say) n x p.
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
#' @param w Censoring weights. This should be a vector of length n. 
#'   Default is NULL, in which case censoring weights are computed using the 
#'   Kaplan-Meier model. 
#' @param active.fun A function which takes the output of train.fun, and reports
#'   which features are active for each fitted model contained in this output.
#'   Its only input argument should be out: output produced by train.fun.
#' @param vars A list specifying the variables (indices between 1 and p)
#'   for which variable importance should be investigated. Alternatively, if
#'   set equal to 0, the default, then all variables are investigated. If both 
#'   vars and active.fun are passed, the latter takes priority and the former
#'   is ignored.
#' @param alpha Miscoverage level for the confidence intervals, i.e., intervals
#'   with coverage 1-alpha are formed. Default for alpha is 0.1.
#' @param bonf.correct Should a Bonferroni correction be applied to the p-values
#'   and confidence intervals? Default is FALSE.
#' @param split Indices that define the data-split to be used (i.e., the indices
#'   define the first half of the data-split, on which the model is trained).
#'   Default is NULL, in which case the split is chosen randomly.
#' @param seed Integer to be passed to set.seed before defining the random 
#'   data-split to be used. Default is NULL, which effectively sets no seed.
#'   If both split and seed are passed, the former takes priority and the latter
#'   is ignored.
#' @param verbose Should intermediate progress be printed out? Default is FALSE.
#'
#' @return A list with the following components: inf.sign, active, master, 
#'   bonf.correct. The first is a list, containing the results of LOCO inference
#'   with the sign text adapted to right-censored data. This list has one 
#'   element per regression model the training and prediction functions, 
#'   train.fun and predict.fun. The inherent to second returned component active
#'   is a list, with one element per tuning step, that reports which features 
#'   are active (or at least investigated) in  the corresponding fitted model.
#'   The third returned component master collects all active (or investigated) 
#'   features across all regression models, for easy reference. The last
#'   component signals whether a Bonferroni correction has been applied.
#'
#' @details The sign test adapted to right-censored data is asymptotically valid
#'   if the Kaplan-Meier model is used to estimate the censoring survival 
#'   function in the censoring weights, under the assumption of independent 
#'   censoring.
#'
#' @references "Distribution-Free Predictive Inference for Regression" by Lei,
#'   G'Sell, Rinaldo, Tibshirani, Wasserman (2018). See also "A Comprehensive 
#'   Framework for Evaluating Time to Event Predictions using the Restricted 
#'   Mean Survival Time" by Cwiling, Perduca, Bouaziz (2023) for the extension
#'   to right-censored data.
#' @export loco.surv

loco.surv = function(x, t, d, tau, train.fun, predict.fun, w=NULL,
                     active.fun=NULL, vars=0, alpha=0.1, rho=0.5, bonf.correct=FALSE, split=NULL, 
                     seed=NULL, verbose=FALSE) {
  
  # Set up data
  # x = as.matrix(x)
  t = as.numeric(t)
  d = as.numeric(d)
  if(!is.null(w)) w = as.numeric(w)
  n = nrow(x)
  p = ncol(x)
  
  # Check input arguments
  check.args.surv(x=x,t=t,d=d,tau=tau,x0=x,w=w,cens.model="km",
                  alpha=alpha,train.fun=train.fun,predict.fun=predict.fun)
  
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
    i1 = sample(1:n,floor(n*rho))
  }
  i2 = (1:n)[-i1]
  n1 = length(i1)
  n2 = length(i2)
  
  if (verbose) {
    cat(sprintf("%sSplitting data into parts of size %i and %i ...\n",
                txt,n1,n2))
    cat(sprintf("%sInitial training on first part ...\n",txt))
  }
  
  # Censoring weights
  if(is.null(w)){
    if (verbose) cat(sprintf("%sComputing censoring weights ...\n",txt))
    w = ipcw(t,d,x,tau,cens.model="km")
  }
  
  # Train on first part
  out = train.fun(x[i1,,drop=F],t[i1],d[i1],tau)
  
  # Which variables to try dropping?
  if(is.null(active.fun)){
    if (length(vars) == 1 && vars == 0) vars = 1:p
    else {
      if (is.null(vars) || !is.numeric(vars) || min(vars) < 1 || 
          max(vars) > p || vars != round(vars)) {
        stop("vars must be a list of integers between 1 and p")
      }
      vars = unique(vars)
    }
    active = purrr::map(1:length(out), ~ vars)
  }else{
    active = active.fun(out)
  }
  
  master = unique(unlist(active))
  m = length(active)
  J = length(master)
  
  that2 = matrix(predict.fun(out,x[i2,,drop=F],tau),nrow=n2)
  
  if (verbose) {
    cat(sprintf("%sInitial residuals on second part ...\n",txt))
  }
  
  # Get residuals on second
  res = abs(pmin(t[i2],tau) - that2)
  res.drop = vector(mode="list",length=J)
  
  
  # Re-fit after dropping each feature in the master list
  for (j in Seq(1,J)) {
    if (verbose) {
      cat(sprintf(paste("\r%sRe-fitting after dropping feature %i (%i of %i",
                        "considered) ..."),txt,master[j],j,J))
      flush.console()
    }
    
    # Train on the first part, without variable 
    out.j = train.fun(x[i1,-master[j],drop=F],t[i1],d[i1],tau)
    
    # Get predictions on the other, without variable
    pred.j = matrix(predict.fun(out.j,x[i2,-master[j],drop=F],tau),nrow=n2)
    res.drop[[j]] = abs(pmin(t[i2],tau) - pred.j)  
  }
  if (verbose) cat("\n")
  
  # Compute p-values and confidence intervals at each tuning step
  inf.sign = vector(mode="list",length=m)
  
  for (l in 1:m) {
    if (verbose) {
      if (m>1) {
        cat(sprintf(paste("\r%sPerforming LOCO analysis at tuning step %i",
                          "(of %i) ..."),txt,l,m))
        flush.console()
      }
      else cat(sprintf("%sPerforming LOCO analysis ...",txt))
    }
    
    k = length(active[[l]])
    inf.sign[[l]] = matrix(0,k,3)
    for (j in Seq(1,k)) {
      j.master = which(master == active[[l]][j])
      z = res.drop[[j.master]][,l] - res[,l]
      inf.sign[[l]][j,] = my.surv.sign.test(t[i2],d[i2],w[i2],z,
                                            tau,alpha,k,bonf.correct)
    }
    
    rownames(inf.sign[[l]]) = active[[l]]
    colnames(inf.sign[[l]]) = c("P-value", "LowConfPt", "UpConfPt")
  }
  if (verbose) cat("\n")
  
  out = list(inf.sign=inf.sign, 
             active=active, master=master, bonf.correct=bonf.correct)
  class(out) = "loco.surv"
  return(out)
}


#' Print function for loco.surv object.
#'
#' @param x The loco.surv object.
#' @param digits Number of digits to display. Default is 3. 
#' @param ... Other arguments (currently not used).
#'
#' @export 

print.loco.surv = function(x, digits=3, ...) {
  
  cat(paste("\nDisplaying results of LOCO analysis. Notes:\n",
            "\n- P-values are from a one-sided test of the target parameter",
            "\n  (mean or median excess test error) being equal to zero versus",
            "\n  greater than zero.",
            "\n- Confidence intervals are from inverting the two-sided version",
            "\n  of this test.",
            "\n- Bonferonni correction:",
            ifelse(x$bonf.correct, "yes.", "no."), "\n"))
  
  m = length(x$inf.sign)
  for (i in 1:m) {
    if (m>1) cat(sprintf("\n------- Tuning step %i -------\n",i))
    
    cat("\nSign test results:\n\n")
    tab = round(x$inf.sign[[i]],digits=digits)
    code = get.signif.code(x$inf.sign[[i]][,1])
    tab = cbind(rownames(tab),tab,code)
    colnames(tab)[1] = "Var"
    colnames(tab)[5] = ""
    rownames(tab) = rep("",nrow(tab))
    print(tab,quote=FALSE)
    
    cat("\nSignificance codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  }
}



# Median inference: one-sided p-value but a two-sided confidence interval, 
# for right-censored data
my.surv.sign.test = function(t, d, w, z, tau, alpha, k, bonf.correct=TRUE){
  n = length(t)
  phi = (z>0)*1 + (z==0)*rbinom(n,1,0.5) 
  w[t > tau] = 0
  
  fit = survfit(Surv(t, d) ~ 1, data=data.frame(t,d))
  S.tau = fit$surv[max(which(fit$time<=tau))]
  
  df.sort = arrange(data.frame(t, d, w, phi),t) %>%
    mutate(Y.Ti = (n-1:n+1)/n,
           phi.w = phi*w,
           phi.bar = rev(cumsum(rev(phi.w)))/n,
           squared.diff = (phi.w - d*phi.bar/Y.Ti - (t<=tau)*d*S.tau/(1-S.tau)*mean(phi*w)/Y.Ti)**2
    )
  
  V = mean(df.sort$squared.diff)/(1-S.tau)**2
  
  tn2 = sqrt(n/V)*(mean(df.sort$phi.w)/(1-S.tau) - 0.5)
  pval = 1-pnorm(tn2)
  
  # Apply Bonferroni correction for k tests
  if (bonf.correct) {
    pval = min(k*pval,1)
    alpha = alpha/k
  }
  
  left = mean(df.sort$phi.w)/(1-S.tau) - sqrt(V/n)*qnorm(1-alpha/2)
  right = mean(df.sort$phi.w)/(1-S.tau) + sqrt(V/n)*qnorm(1-alpha/2)
  
  return(c(pval,left,right))
}