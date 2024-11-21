#' Evaluating predictions of a Restricted Mean Survival Time estimator
#' 
#' Evaluation of a learning algorithm for the restricted mean survival time
#' trained on right-censored data. Estimation of its mean squared error, 
#' computation of prediction intervals, and evaluation of (local and global) 
#' variable importance.
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
#'   training for several (say m) learning models.
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
#' @param CV Boolean indicating whether or not to perform cross-validation for
#'   the estimation of the mean squared error. Default is TRUE.
#' @param n.folds The number of folds for the cross-validated estimation of the
#'   mean squared error. The default number is 10 folds.
#' @param prop Proportion of data to put in the test set in case cross-
#'   validation is not performed. If prop=0 then the whole data set is used 
#'   both as training set and test set.
#' @param active.fun A function which takes the output of train.fun, and reports
#'   which features are active for each fitted model contained in this output.
#'   Its only input argument should be out: output produced by train.fun.
#' @param alpha Miscoverage level for the prediction intervals, i.e., intervals
#'   with coverage 1-alpha are formed. Default for alpha is 0.1.
#' @param varsL A list specifying the variables (indices between 1 and p)
#'   for which local variable importance should be investigated. Alternatively,
#'   if set equal to 0, the default, then all variables are investigated.
#' @param varsG A list specifying the variables (indices between 1 and p)
#'   for which global variable importance should be investigated. Alternatively,
#'   if set equal to 0, the default, then all variables are investigated.
#' @param bonf.correct Should a Bonferroni correction be applied to the p-values
#'   and confidence intervals? Default is FALSE.
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
#'   to conformal.pred.roo.surv (whether already run, or to be run inter-
#'   nally, see below), loco.roo.surv and loco.surv.
#'   Default is NULL, in which case the split is chosen randomly.
#' @param rho Proportion of the data for the data-split in case it is chosen
#'   randomly.
#' @param seed Integer to be passed to set.seed. Default is NULL, which 
#'   effectively sets no seed.
#' @param out.roo.surv Output from running conformal.pred.roo.surv on the given  
#'   data set; default is NULL, which means that conformal.pred.roo.surv will be  
#'   run internally in order to compute in-sample prediction intervals. If  
#'   is NULL, then the data-split common to conformal.pred.roo.surv and this 
#'   out.roo.surv function will be determined by the split or seed arguments. If  
#'   non-NULL, then the common data-split will be determined by that recorded in  
#'   out.roo.surv, and the split argument will be (silently) ignored.
#' @param verbose Should intermediate progress be printed out? Default is FALSE.
#' @param error Boolean indicating whether or not to estimate the mean squared 
#' error. Default is TRUE.
#' @param roo Boolean indicating whether or not to compute in-sample rank-one-
#' out prediction intervals. Default is TRUE.
#' @param vimpL Boolean indicating whether or not to compute local variable
#' importance. Default is TRUE.
#' @param vimpG Boolean indicating whether or not to compute global variable
#' importance. Default is TRUE.
#'
#' @return A list with the following components: x, split, mse, out.roo.surv, 
#'   out.loco.roo.surv, out.loco.surv. The matrix of variables x is returned in
#'   x. The indices used for the half of the data-split are returned in split. 
#'   mse is a matrix of dimension n.folds x m containing the estimation of the 
#'   mean squared error using censoring weights on each fold and for each 
#'   learning model. The last three components are the outputs of respectively
#'   roo.surv, loco.roo.surv and loco.surv.
#'
#' @details Analysis framework of one or more restricted mean survival time
#'   estimation model(s). Details are given in the descriptions of functions
#'   wrss, conformal.pred.roo.surv, loco.roo.surv and loco.surv. The bias of 
#'   computing residuals on observed censored times is corrected by using 
#'   censoring weights. Convergence guarantees for global variable importance 
#'   currently hold only if the Kaplan-Meier estimator is used to compute 
#'   censoring weights.
#'   
#' @seealso \code{\link{wrss}}, \code{\link{conformal.pred.roo.surv}}, 
#'   \code{\link{loco.roo.surv}}, \code{\link{loco.surv}}
#'
#' @references See "Distribution-Free Predictive Inference for Regression" by 
#'   Lei, G'Sell, Rinaldo, Tibshirani, Wasserman (2018) as a reference for 
#'   conformal prediction. See "A Comprehensive Framework for Evaluating Time to 
#'   Event Predictions using the Restricted Mean Survival Time" by Cwiling, 
#'   Perduca, Bouaziz (2023) for the extension to right-censored data.
#' @export rmst.pred

rmst.pred = function(x, t, d, tau,  train.fun, predict.fun, w=NULL, cens.model="km", 
                     CV=T, n.folds=10, prop=0.1, active.fun=NULL, alpha=0.1,  
                     varsL=0, varsG=0, bonf.correct=FALSE, mad.train.fun=NULL, 
                     mad.predict.fun=NULL, split=NULL, rho=0.5, seed=NULL, 
                     out.roo.surv=NULL,  verbose=FALSE, error=T,roo=T,vimpL=T,vimpG=T) {
  
  # Set up data
  # x = as.matrix(x)
  t = as.numeric(t)
  d = as.numeric(d)
  if(!is.null(w)) w = as.numeric(w)
  n = nrow(x)
  p = ncol(x)
  
  # Check input arguments
  check.args.surv(x=x,t=t,d=d,tau=tau,x0=x,w=w,cens.model=cens.model,
                  alpha=alpha,train.fun=train.fun,predict.fun=predict.fun,
                  mad.train.fun=mad.train.fun,mad.predict.fun=mad.predict.fun)
  
  # List of elements to return
  out = list(m=NULL,x=x,split=NULL,mse=NULL,out.roo.surv=NULL,
             out.loco.roo.surv=NULL, out.loco.surv=NULL)
  
  # Users may pass in a string for the verbose argument
  if (verbose == TRUE) txt = ""
  if (verbose != TRUE && verbose != FALSE) {
    txt = verbose
    verbose = TRUE
  }
  
  # Seed
  if (!is.null(seed)) set.seed(seed)
  
  # Censoring weights
  if(is.null(w)){
    if (verbose) cat(sprintf("%sComputing censoring weights ...\n",txt))
    w = ipcw(t,d,x,tau,cens.model)
  }
  
  if(error){
    # Estimate the mean squared error
    if (verbose) {
      cat(sprintf(paste("%sComputing mean squared error ...\n"),txt))
    }
    mse = wrss(x,t,d,tau,train.fun,predict.fun,w=w,cens.model=cens.model,
               CV=CV,n.folds=n.folds,prop=prop)
    out$mse = mse
    out$m = ifelse(is.null(nrow(mse)),length(mse),ncol(mse))
  }
  
  
  # Compute in-sample rank-one-out intervals ourselves
  if (is.null(out.roo.surv)) {
    # If the user passed indices for the split, use them, 
    # otherwise make a random split
    if (is.null(split)) split = sample(1:n,floor(n*rho))
    
    if(roo){
      if (verbose) {
        cat(sprintf(paste("%sComputing in-sample rank-one-out conformal",
                          "prediction intervals ...\n"),txt))
      }
      
      out.roo.surv = conformal.pred.roo.surv(x,t,d,tau,train.fun,predict.fun,
                                             alpha=alpha,mad.train.fun=mad.train.fun,mad.predict.fun=mad.predict.fun, 
                                             split=split,w=w, cens.model=cens.model,verbose=F)
      out$out.roo.surv=out.roo.surv
      out$m = ncol(out.roo.surv$pred)
    }
  }
  # These have already been computed
  else {
    split = out.roo.surv$split
  }
  out$split = split
  
  
  if(vimpL){
    # Local variable importance
    if (verbose) {
      cat(sprintf(paste("%sComputing local variable importance ...\n"),txt))
    }
    out.loco.roo.surv = loco.roo.surv(x, t, d, tau, train.fun, predict.fun, 
                                      vars=varsL, alpha=alpha, out.roo.surv=out.roo.surv, verbose=F)
    out$out.loco.roo.surv=out.loco.roo.surv
    out$m = ncol(out.loco.roo.surv$out.roo.surv$pred)
  }
  
  if(vimpG){
    # Global variable importance
    if (verbose) {
      cat(sprintf(paste("%sComputing global variable importance ...\n"),txt))
    }
    out.loco.surv = loco.surv(x, t, d, tau, train.fun, predict.fun, 
                              active.fun=active.fun, vars=varsG, alpha=alpha, bonf.correct=bonf.correct, 
                              split=split, verbose=F)
    out$out.loco.surv=out.loco.surv
    out$m = length(out.loco.surv$active)
  }
  
  
  class(out) = "rmst.pred"
  return(out)
}


#' Print function for rmst.pred object.
#'
#' @param x The rmst.pred object.
#' @param elements Selection among "mse", "intervals", "vimp", or "all",   
#'   describing which elements of the analysis to show. If "mse", then the  
#'   results on the estimation of the mean squared error are displayed; if 
#'   "intervals", then information on the length of the prediction intervals are 
#'   displayed; if  "vimp", then results of global variable importance are 
#'   displayed. If "all", then the results from all elements are displayed. 
#'   Default is "all".
#' @param digits Number of digits to display. Default is 3.
#' @param model.names Names given to the different learning models given to
#'   train.fun. Default is NULL.
#' @param var.names Names given to the different variables used in the 
#' regression. Default is NULL.   
#' @param ... Other arguments (currently not used).
#'
#' @export 

print.rmst.pred = function(x,elements=c("all"),digits=3,
                           model.names=NULL, var.names=NULL, ...) {
  
  n.folds = ifelse(is.null(nrow(x$mse)),1,nrow(x$mse))
  if(is.null(x$m)) stop("Nothing to print")
  m = x$m
  if(is.null(model.names)){model.names = paste("Model",1:m)}
  if(!is.character(model.names) || length(model.names)!=m) 
    stop("model.names must be a vector of strings of size m")
  if(!("mse" %in% elements) && !("intervals" %in% elements) && 
     !("vimp" %in% elements) && !("all" %in% elements))
    stop("elements must contain at least mse, intervals, vimp or all")
  
  if(("mse" %in% elements || "all" %in% elements) && !is.null(x$mse)){
    if(n.folds==1){
      cat("\n*Estimation of the mean squared error.\n")
      mse.disp = paste(round(x$mse,digits=digits))
    }else{
      cat(paste("\n*Estimation of the mean squared error on",nrow(x$mse),"folds:",
                "\n  Mean and quantiles are computed and displayed as follows:",
                "\n  Mean [0.1-quantile, 0.9-quantile]\n"))
      mse.disp = paste(round(apply(x$mse,2,mean),digits=digits),"[",
                       round(apply(x$mse,2,function(s)quantile(s,0.1)),digits=digits),",",
                       round(apply(x$mse,2,function(s)quantile(s,0.9)),digits=digits),"]")
    }
  }
  
  
  if(("intervals" %in% elements || "all" %in% elements) && !is.null(x$out.roo.surv)){
    d = x$out.roo.surv$up - x$out.roo.surv$lo
    
    if(sum(apply(d,2,function(s)length(unique(s)))) == m){
      cat("\n*Length of the prediction intervals (centered on the predictions).\n")
      d.disp = paste(round(apply(d,2,mean),digits=digits))
    }else{
      cat("\n*Mean and quantiles of the length of the prediction intervals",
          "\n  (centered on the predictions) displayed as follows:",
          "\n  Mean [0.1-quantile, 0.9-quantile]\n")
      
      d.disp = paste(round(apply(d,2,mean),digits=digits),"[",
                     round(apply(d,2,function(s)quantile(s,0.1)),digits=digits),",",
                     round(apply(d,2,function(s)quantile(s,0.9)),digits=digits),"]")
    }
  }
  
  
  if(("vimp" %in% elements || "all" %in% elements) && !is.null(x$out.loco.surv)){
    cat(paste("\n*Displaying results of LOCO analysis. Notes:\n",
              "\n - P-values are from a one-sided test of the target parameter",
              "\n  (mean or median excess test error) being equal to zero versus",
              "\n  greater than zero.",
              "\n - Confidence intervals are from inverting the two-sided version",
              "\n  of this test.",
              "\n - Bonferonni correction:",
              ifelse(x$out.loco.surv$bonf.correct, "yes.", "no."), "\n"))
  }
  
  
  
  for (i in 1:m) {
    if (m>1) cat(sprintf("\n------- %s -------\n",model.names[i]))
    
    if(("mse" %in% elements || "all" %in% elements) && !is.null(x$mse)){
      cat("\nMean squared error estimation:\n")
      cat(sprintf("\n %s \n", mse.disp[i]))
    }
    
    if(("length" %in% elements || "all" %in% elements) && !is.null(x$out.roo.surv)){
      cat("\nLength of the prediction interval around the prediction:\n")
      cat(sprintf("\n %s \n", d.disp[i]))
    }
    
    if(("vimp" %in% elements || "all" %in% elements) && !is.null(x$out.loco.surv)){
      cat("\nSign test results:\n\n")
      tab = round(x$out.loco.surv$inf.sign[[i]],digits=digits)
      if(!is.null(var.names)) rownames(tab) = var.names
      code = get.signif.code(x$out.loco.surv$inf.sign[[i]][,1])
      tab = cbind(rownames(tab),tab,code)
      colnames(tab)[1] = "Var"
      colnames(tab)[5] = ""
      rownames(tab) = rep("",nrow(tab))
      print(tab,quote=FALSE)
      
      cat("\nSignificance codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
    }
  }
}


#' Plot function for rmst.pred object.
#'
#' @param x The rmst.pred object.
#' @param elements Selection among "mse", "vimpL", "vimpG" or "all"describing 
#'   which elements of the analysis to show. If "mse" then the results on the 
#'   estimation of the mean squared error are displayed; if "vimpL" then results
#'   of local variable importance are displayed; if "vimpG" then results of 
#'   global variable importance are displayed. If "all" then the results from 
#'   all elements are displayed. Default is "all". 
#' @param model.names Names given to the different learning models given to
#'   train.fun. Default is NULL.
#' @param ... Other arguments (currently not used).
#'
#' @export 

plot.rmst.pred = function(x,elements=c("all"), model.names=NULL, ...) {
  
  n.folds = ifelse(is.null(nrow(x$mse)),1,nrow(x$mse))
  if(is.null(x$m)) stop("Nothing to plot")
  m = x$m
  if(is.null(model.names)){model.names = paste("Model",1:m)}
  if(!is.character(model.names) || length(model.names)!=m) 
    stop("model.names must be a vector of strings of size m")
  if(!("mse" %in% elements) && 
     !("vimpL" %in% elements) && !("vimpG" %in% elements) 
     && !("all" %in% elements))
    stop("elements must contain at least mse, vimpL, vimpG or all")
  
  if(("mse" %in% elements || "all" %in% elements) && !is.null(x$mse)){
    V = nrow(x$mse)
    tab = data.frame(rep=1:V,x$mse)
    colnames(tab) = c("rep",model.names)
    tab = reshape::melt(tab,id="rep")
    colnames(tab) = c("rep","model","mse")
    plot = tab %>% 
      ggplot(aes(x = as.factor(model), y = mse)) +
      geom_boxplot(lwd=0.3,outlier.size = 0.3, fatten=1) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) +
      labs(x = 'Learning model', y = 'Estimated mean squared error')
    print(plot)
  }
  
  if(("vimpL" %in% elements || "all" %in% elements) && !is.null(x$out.loco.roo.surv)){
    varsL = x$out.loco.roo.surv$vars
    for (i in 1:m){
      for (j in 1:length(varsL)){
        plot(x$x[,varsL[j]],x$x[,varsL[j]],ylim=range(c(x$out.loco.roo.surv$lo,
                                                        x$out.loco.roo.surv$up)),xlab="Location",ylab="Interval",
             main=paste(model.names[i], "- Var.",varsL[j]),col=NA)
        cols = ifelse(x$out.loco.roo.surv$lo[,j,i] <= 0, 1, 3)
        segments(x$x[,varsL[j]],x$out.loco.roo.surv$lo[,j,i],
                 x$x[,varsL[j]],x$out.loco.roo.surv$up[,j,i],
                 col=cols,lwd=1)
        abline(h=0, lty=2, lwd=2, col=2)
      }
    }
  }
  
  if(("vimpG" %in% elements || "all" %in% elements) && !is.null(x$out.loco.surv)){
    ylim = range(unlist(map(x$out.loco.surv$inf.sign,~range(.[,2:3]))))
    for (i in 1:m){
      p = length(x$out.loco.surv$active[[i]])
      plot(c(), c(), xlim=c(1,p), ylim=ylim, xaxt="n",
           main=model.names[i],
           xlab="Variable", ylab="Confidence interval")
      axis(side=1, at=1:p, labels=FALSE)
      for (j in 1:p) {
        axis(side=1, at=j, labels=x$out.loco.surv$active[[i]][j], cex.axis=0.75,
             line=0.5*j%%2)
      }
      abline(h=0.5)
      segments(1:p, x$out.loco.surv$inf.sign[[i]][,2],
               1:p, x$out.loco.surv$inf.sign[[i]][,3],
               col="red", lwd=2)
      points(1:p,(x$out.loco.surv$inf.sign[[i]][,3]+x$out.loco.surv$inf.sign[[i]][,2])/2,
             col="red",pch=23)
    }
  }
  
}

#' Mean Squared Error of a Restricted Mean Survival Time estimator
#' 
#' Estimation of the mean squared error of a learning algorithm for the 
#' restricted mean survival time trained on right-censored data. 
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
#'   training for several (say m) learning models.
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
#' @param CV Boolean indicating whether or not to perform cross-validation.
#'   Default is TRUE.
#' @param n.folds The number of folds for the cross-validated estimation of the
#'   mean squared error. The default number is 10 folds.
#' @param split Indices that define the data-split to be used (i.e., the indices
#'   define the training set) in case cross-validation is not performed. 
#'   Default is NULL, in which case the split is chosen randomly.
#' @param prop Proportion of data to put in the test set in case cross-
#'   validation is not performed. If prop=0 then the whole data set is used 
#'   both as training set and test set.
#' @param seed Integer to be passed to set.seed. Default is NULL, which 
#'   effectively sets no seed.
#'
#' @return A matrix of dimension n.folds (1 if cross-validation is not 
#'   performed) x m containing the estimation of the mean squared error using 
#'   censoring weights on each fold and for each learning model. 
#'
#' @details 
#' The estimation of the mean squared error with the WRSS estimator
#' can be performed in one of the following ways: 
#' 
#' - By cross-validation on a defined number of folds;
#'
#' - By dividing the data into a training set and a test set, on which to train
#'   the learning algorithm and compute the estimator of the mean squared error,
#'   respectively;
#'
#' - By using all data both as a training set and a test set.
#'
#' The default method is the first one.
#'
#' @references See "A Comprehensive Framework for Evaluating Time to 
#'   Event Predictions using the Restricted Mean Survival Time" by Cwiling, 
#'   Perduca, Bouaziz (2023).
#' @export wrss

wrss = function(x,t,d,tau,train.fun,predict.fun,w=NULL,cens.model="km",
                CV=T,n.folds=10,split=NULL,prop=0.1,seed=NULL){
  check.args.surv(x=x,t=t,d=d,tau=tau,x0=x,w=w,cens.model=cens.model,
                  alpha=0.1,train.fun=train.fun,predict.fun=predict.fun)
  check.bool(CV)
  
  n = length(t)
  if (is.null(n.folds) || !is.numeric(n.folds) || length(n.folds) > 1 ||  
      n.folds < 2 || n.folds > n || n.folds != round(n.folds)) {
    stop("n.folds must be an integer between 2 and n")
  }
  
  check.num.01(prop)
  
  if (!is.null(seed)) set.seed(seed)
  
  if(is.null(w)){w = ipcw(t,d,x,tau,cens.model)}
  
  if(CV){
    # Cross-validation
    folds = sample(cut(1:n,breaks=n.folds,labels=FALSE))
    folds = purrr::map(1:n.folds, ~ which(folds==.,arr.ind=TRUE))
    trainings = purrr::map(folds, ~ train.fun(x[-.,],t[-.],d[-.],tau))
    fits = purrr::map2(folds,trainings, ~ matrix(predict.fun(.y,x[.x,],tau),nrow=length(.x)))
    mse = purrr::map2(folds,fits, ~ apply(w[.x]*(pmin(t[.x],tau) - .y)**2,2,sum)/length(w[.x]))
    mse = matrix(unlist(mse),nrow=n.folds,byrow=T)
  }else{
    if(prop==0){
      # Train and fit on full data set
      out.all = train.fun(x,t,d,tau) 
      fit.all = matrix(predict.fun(out.all,x,tau),nrow=n)
      mse = apply(w*(pmin(t,tau) - fit.all)**2,2,sum)/n
    }else{
      # If the user passed indices for the split, use them
      if (!is.null(split)) i1 = split
      # Otherwise make a random split
      else {
        i1 = sample(1:n,floor(n*(1-prop)))
      }
      i2 = (1:n)[-i1]
      # Train on the first data set and fit on the second
      out.all = train.fun(x[i1,,drop=F],t[i1],d[i1],tau) 
      fit.all = matrix(predict.fun(out.all,x[i2,,drop=F],tau),nrow=length(i2))
      mse = apply(w[i2]*(pmin(t[i2],tau) - fit.all)**2,2,sum)/length(i2)
    }
  }
  return(mse)
}

