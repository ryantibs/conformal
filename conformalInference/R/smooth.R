#' Smoothing spline training and prediction functions, and special leave-one-out
#'   fitting function.
#'
#' Construct the training and prediction functions for smoothing splines, at a
#'   particular tuning parameter value or df value (either given or chosen by
#'   cross-validation), and construct also a special function for computing
#'   leave-one-out fitted values. Based on the \code{\link{smooth.spline}}
#'   function.
#'
#' @param df Desired effective degrees of freedom of the smoothing spline fit.
#'   If NULL, then the complexity of the fit controlled either the spar or cv
#'   arguments (see below).
#' @param spar Parameter that scales roughly as log(lambda), where lambda is
#'   the parameter multiplying the roughness penalty in the smoothing spline
#'   criterion. See the smooth.spline help file for details. If NULL, then the
#'   complexity of the fit is chosen by the cv argument (see next).
#' @param cv Either TRUE, indicating that leave-one-out cross-validation should
#'   be used to choose the smoothing parameter lambda, or FALSE, indicating that
#'   generalized cross-validation should be used to choose lamdbda. Note that
#'   this argument is only in effect when both df and spar are NULL.
#' @param tol.fact Factor multiplying the interquartile range of the x values,
#'   IQR(x), to determine a threshold for uniqueness: the x values are binned
#'   into bins of size tol.fact * IQR(x), and values which fall into the same
#'   bin are regarded as the same. Default is 1e-6.
#'
#' @details Note that the constructed function special.fun only serve as an
#'   approximation to the exact leave-one-out fits when the smooth.spline
#'   function merges nearby adjacent x points (which happens when adjacent x
#'   points are closer than tol), since in this case the smooth.spline function
#'   only returns the leverage scores (diagonal elements of the smoother matrix)
#'   on the reduced (merged) data set.  Another important note is that the
#'   constructed special.fun will be formally invalid when cross-validation (or
#'   generalized cross-validation) is used to pick the smooth parameter. In
#'   both of these cases, special.fun does \emph{not} throw a warning, and it is
#'   left up to the user to consider using this function wisely.
#'
#' @return A list with four components: train.fun, predict.fun, special.fun,
#'   active.fun. The last function is designed to take the output of train.fun,
#'   and reports which features are active for each fitted model contained in
#'   this output. Trivially, here, there is only one feature and it is always
#'   active.
#'
#' @author Ryan Tibshirani
#' @example examples/ex.smooth.funs.R
#' @export smooth.spline.funs

smooth.spline.funs = function(df=NULL, spar=NULL, cv=TRUE, tol.fact=1e-6) {
  if (!is.null(df)) check.pos.int(df)
  if (!is.null(spar)) check.pos.num(spar)
  check.bool(cv)
  check.pos.num(tol.fact)
  
  if (!is.null(df)) {
    train.fun = function(x,y) {
      return(suppressWarnings(smooth.spline(x,y,df=df,spar=spar,
                                            tol=tol.fact*IQR(x))))
    }
  }
  else {
    train.fun = function(x,y) {
      return(suppressWarnings(smooth.spline(x,y,spar=spar,cv=cv,
                                            tol=tol.fact*IQR(x))))
    }
  }

  predict.fun = function(out,newx) {
    return(predict(out,newx)$y)
  }

  special.fun = function(x,y,out) {
    n = nrow(x)
    if (length(out$lev) == n) s = out$lev
    else s = approx(out$x,out$lev,xout=x)$y
    return((y - predict(out,x)$y) / (1-s))
  }

  active.fun = function(out) {
    return(list(1))
  }
  
  return(list(train.fun=train.fun, predict.fun=predict.fun,
              special.fun=special.fun, active.fun=active.fun))
}
