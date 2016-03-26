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
#' Some key references (in reverse chronological order):
#' - "Distribution-Free Predictive Inference for Regression" by Max G'Sell,
#'   Jing Lei, Alessandro Rinaldo, Ryan Tibshirani, and Larry Wasserman,
#'   http://arxiv.org/pdf/xxxx.pdf, 2016.
#' - "Classification with Confidence" by Jing Lei (2014).
#' - "Distribution-Free Prediction Bands for Non-parametric Regression" by
#'   Jing Lei and Larry Wasserman (2014). 
#' - "A Conformal Prediction Approach to Explore Functional Data" by Jing Lei,
#'   Alessandro Rinaldo, and Larry Wasserman (2013).
#' - "Distribution Free Prediction Sets" by Jing Lei, James Robins, and Larry 
#'   Wasserman (2013).
#' - "On-line Predictive Linear Regression" by Vladimir Vovk, Ilia Nouretdinov,
#'   and Alex Gammerman (2009).
#' - "Algorithmic Learning in a Random World" by Vladimir Vovk, Alex Gammerman,
#'   and Glenn Shafer (2005). 
"_PACKAGE"
#> [1] "_PACKAGE"