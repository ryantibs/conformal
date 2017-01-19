#' Tools for conformal inference in regression
#'
#' A collection of tools for distribution-free inference for regression problems
#'   using the theory of conformal prediction. 
#'
#' Conformal inference is a framework for converting any pre-chosen estimator of
#'   the regression function into prediction intervals with finite-sample
#'   validity, under essentially no assumptions on the data-generating process
#'   (aside from the the assumption of i.i.d. observations). The main functions
#'   in this package for computing such prediction intervals are
#'   \code{\link{conformal.pred}} (shorter intervals, but more expensive) and
#'   \code{\link{conformal.pred.split}} (wider intervals, but much more
#'   efficient).
#'
#' Main reference: "Distribution-Free Predictive Inference for Regression" by
#'   Jing Lei, Max G'Sell, Alessandro Rinaldo, Ryan Tibshirani, and Larry
#'   Wasserman, http://arxiv.org/pdf/1604.04173.pdf, 2016.
#'
#' Related references (in reverse chronological order):
#' \itemize{
#' \item "Classification with Confidence" by Jing Lei, Biometrika, 101(4), 
#'   755-769, 2014. 
#' \item "Distribution-Free Prediction Bands for Non-parametric Regression" by
#'   Jing Lei and Larry Wasserman, Journal of the Royal Statistical Society: 
#'   Series B, 76(1), 71-96, 2014. 
#' \item "A Conformal Prediction Approach to Explore Functional Data" by Jing 
#'   Lei, Alessandro Rinaldo, and Larry Wasserman, Annals of Mathematics and 
#'   Artificial Intelligence, 74(4), 29-43, 2013.  
#' \item "Distribution Free Prediction Sets" by Jing Lei, James Robins, and 
#'   Larry  Wasserman, Journal of the American Statistical Association, 
#'   108(501), 278-287, 2013.  
#' \item "On-line Predictive Linear Regression" by Vladimir Vovk, Ilia 
#'   Nouretdinov, and Alex Gammerman, Annals of Statistics, 37(3), 1566-1590,
#'   2009.
#' \item "Algorithmic Learning in a Random World" by Vladimir Vovk, Alex 
#'   Gammerman, and Glenn Shafer, Springer, 2005.
#' }
"_PACKAGE"
#> [1] "_PACKAGE"
