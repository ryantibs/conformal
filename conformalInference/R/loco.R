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
#'   active, master. The first three are lists, containing the results of LOCO
#'   inference with the Z-test, sign text, and Wilcoxon signed rank test,
#'   respectively. More
#'   details on these tests are given below. These lists have one element per
#'   tuning step inherent to the training and prediction functions, train.fun
#'   and predict.fun. The fourth returned component active is a list, with one
#'   element per tuning step, that reports which features are active in the
#'   corresponding fitted model. The last returned component master collects
#'   all active features across all tuning steps, for easy reference.
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
#'   A few other important notes: p-values here are from a one-sided test of the
#'   target parameter (mean or median excess test error) being equal to zero
#'   versus greater than zero. Confidence intervals are from inverting the
#'   two-sided version of this test. Furthermore, all p-values and confidence
#'   intervals have been Bonferroni-corrected for multiplicity.
#'
#' @author Ryan Tibshirani, Larry Wasserman
#' @references "Distribution-Free Predictive Inference for Regression" by 
#'   Jing Lei, Max G'Sell, Alessandro Rinaldo, Ryan Tibshirani, and Larry
#'   Wasserman, https://arxiv.org/pdf/1604.04173.pdf, 2016.
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
  for (j in 1:J) {
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
                          "(of %i) ...",txt,l,m)))
        flush.console()
      }
      else cat(sprintf("%sPerforming LOCO analysis ...",txt))
    }

    k = length(active[[l]])
    inf.z[[l]] = inf.sign[[l]] = inf.wilcox[[l]] = matrix(0,k,3)
    for (j in 1:k) {
      j.master = which(master == active[[l]][j])
      z = res.drop[[j.master]][,l] - res[,l]
      inf.z[[l]][j,] = my.z.test(z,alpha,k)
      inf.sign[[l]][j,] = my.sign.test(z,alpha,k)
      inf.wilcox[[l]][j,] = my.wilcox.test(z,alpha,k)
      
      rownames(inf.z[[l]]) = rownames(inf.sign[[l]]) =
        rownames(inf.wilcox[[l]]) = active[[l]]
      colnames(inf.z[[l]]) = colnames(inf.sign[[l]]) =
        colnames(inf.wilcox[[l]]) = c("P-value", "LowConfPt", "UpConfPt")
    }
  }
  if (verbose) cat("\n")

  out = list(inf.z=inf.z,inf.sign=inf.sign,inf.wilcox=inf.wilcox,
    active=active, master=master)
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
            "\n- All p-values and confidence intervals have been Bonferroni",
            "\n  corrected for multiplicity.\n"))
  
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
my.z.test = function(z, alpha, k){
  n = length(z)
  s = sd(z)
  m = mean(z)
  pval = 1-pnorm(m/s*sqrt(n))

  # Apply Bonferroni correction for k tests
  pval = min(k*pval,1)
  alpha = alpha/k
  
  q = qnorm(1-alpha/2)
  left  = m - q*s/sqrt(n)
  right = m + q*s/sqrt(n)
  return(c(pval,left,right))
}

# Median inference: one-sided p-value but a two-sided confidence interval
my.sign.test = function(z, alpha, k){
  n = length(z)
  s = sum(z>0)
  pval = 1-pbinom(s-1,n,0.5)

  # Apply Bonferroni correction for k tests
  pval = min(k*pval,1)
  alpha = alpha/k
  
  j = n - qbinom(1-alpha/2,n,0.5) # Want P(B <= j) <= alpha/2
  zz  = sort(z)
  left  = zz[j]
  right = zz[n-j+1]
  return(c(pval,left,right))
}

# Median inference: one-sided p-value but a two-sided confidence interval
my.wilcox.test = function(z, alpha, k){
  pval = wilcox.test(z, alternative="greater")$p.value

  # Apply Bonferroni correction for k tests
  pval = min(k*pval,1)
  alpha = alpha/k
  
  out = wilcox.test(z, conf.int=TRUE, conf.level=1-alpha)
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


