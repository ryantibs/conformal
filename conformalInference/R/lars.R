#' Forward stepwise, least angle regression, and lasso training and prediction
#'   functions. 
#'
#' Construct training and prediction functions for forward stepwise, least angle
#'   regression net, or the lasso path, based on the \code{\link{lars}} package,
#'   over a given number of steps.
#'
#' @param type One of "stepwise", "lar", or "lasso", indicating which path
#'   algorithm should be run. Default is "stepwise".
#' @param normalize,intercept Should the features be normalized (to have
#'   unit L2 norm), and should an intercept be included? Default for both is
#'   TRUE.
#' @param max.steps Number of steps down the path to be taken before quitting.
#'   Default is 20. Prediction is done over steps 1 through max.steps, unless
#    cv is TRUE (see next).
#' @param cv Should 10-fold cross-validation be used? If TRUE, then prediction
#'   is done with either the (usual) min rule, or the 1se rule. See the cv.rule
#'   argument, below. If FALSE (the default), then prediction is done over all 
#'   path steps taken.
#' @param cv.rule If the cv argument is TRUE, then cv.rule determines which rule
#'   should be used for the predict function, either "min" (the usual rule) or
#'   "1se" (the one-standard-error rule). Default is "min".
#' @param use.Gram Should the Gram matrix (cross product of feature matrix)
#'   be stored? Default is TRUE. But storing the Gram matrix may be costly when
#'   p (number of colums of feature matrix) is large, say, large compared to n
#'   (number rows of feature matrix), so in this case, you may want to choose
#'   FALSE.
#' 
#' @return A list with three components: train.fun, predict.fun, active.fun.
#'   The third function is designed to take the output of train.fun, and
#'   reports which features are active for each fitted model contained in
#'   this output.
#'
#' @details This function is based on the packages \code{\link{lars}} and
#'   \code{\link{plyr}}. If these packages are not installed, then the function
#'   will abort. 
#'
#' @seealso \code{\link{glmnet.funs}} for elastic net, lasso, and ridge
#'   regression training and prediction functions, defined over a sequence
#'   of lambda values.
#' @author Ryan Tibshirani
#' @example examples/ex.lars.funs.R
#' @export lars.funs

lars.funs = function(type=c("stepwise","lar","lasso"), normalize=TRUE,
  intercept=TRUE, max.steps=20, cv=FALSE, cv.rule=c("min","1se"),
  use.Gram=TRUE) {

  # Check for lars
  if (!require("lars",quietly=TRUE)) {
    stop("Package lars not installed (required here)!") 
  }
  # Check for plyr
  if (!require("plyr",quietly=TRUE)) {
    stop("Package plyr not installed (required here)!")
  }

  # Check arguments
  type = match.arg(type)
  check.bool(normalize)
  check.bool(intercept)
  check.pos.int(max.steps)
  check.bool(cv)
  cv.rule = match.arg(cv.rule)
  check.bool(use.Gram)

  # They want training to be done over a bunch of steps, and prediction at the
  # step number chosen by either the min CV rule or 1se CV rule
  if (cv) {
    train.fun = function(x,y,out=NULL) {
      # Subtlety: trim max.steps if we need to
      max.steps = min(max.steps,nrow(x),ncol(x)-intercept)
      # Fit the lars path object
      out = lars(x,y,type=type,normalize=normalize,intercept=intercept,
        max.steps=max.steps,use.Gram=use.Gram)
      # Fit the lars cv object (inefficient! but we have to)
      out.cv = cv.lars(x,y,index=1:max.steps,mode="step",plot.it=FALSE,
        max.steps=max.steps,normalize=normalize,intercept=intercept,
        use.Gram=use.Gram)
      # Compute k.min and k.1se, the min and 1se choices from CV
      k.min = which.min(out.cv$cv)
      k.1se = min(which(out.cv$cv <= out.cv$cv[k.min] +
        out.cv$cv.error[k.min]))
      # Append k.min and k.1se to the original lars object
      out$k.min = k.min
      out$k.1se = k.1se
      return(out)
    }
    
    predict.fun = function(out,newx) {
      return(predict(out,newx,type="fit",mode="step",
                     s=ifelse(cv.rule=="min",out$k.min,out$k.1se))$fit)
    }

    active.fun = function(out) {
      b = coef(out,mode="step",s=ifelse(cv.rule=="min",out$k.min,out$k.1se))
      return(list(which(b!=0)))
    }
  }

  # They want training to be done over a bunch of steps, and prediction at the
  # same steps
  else {
    train.fun = function(x,y,out=NULL) {
      # Subtlety: trim max.steps if we need to
      max.steps = min(max.steps,nrow(x),ncol(x)-intercept)
      return(lars(x,y,type=type,normalize=normalize,intercept=intercept,
        max.steps=max.steps,use.Gram=use.Gram))
    }

    predict.fun = function(out,newx) {
      return(predict(out,newx,type="fit",mode="step")$fit)
    }

    active.fun = function(out) {
      b = coef(out,mode="step")
      return(alply(b,1,.fun=function(v) which(v!=0)))
    }
  }
  
  return(list(train.fun=train.fun, predict.fun=predict.fun,
              active.fun=active.fun))
}
