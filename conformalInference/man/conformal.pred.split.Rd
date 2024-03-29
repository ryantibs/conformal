% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/split.R
\name{conformal.pred.split}
\alias{conformal.pred.split}
\title{Split conformal prediction intervals.}
\usage{
conformal.pred.split(
  x,
  y,
  x0,
  train.fun,
  predict.fun,
  alpha = 0.1,
  rho = 0.5,
  w = NULL,
  mad.train.fun = NULL,
  mad.predict.fun = NULL,
  split = NULL,
  seed = NULL,
  verbose = FALSE
)
}
\arguments{
\item{x}{Matrix of features, of dimension (say) n x p.}

\item{y}{Vector of responses, of length (say) n.}

\item{x0}{Matrix of features, each row being a point at which we want to
form a prediction interval, of dimension (say) n0 x p.}

\item{train.fun}{A function to perform model training, i.e., to produce an
estimator of E(Y|X), the conditional expectation of the response variable
Y given features X. Its input arguments should be x: matrix of features,
and y: vector of responses.}

\item{predict.fun}{A function to perform prediction for the (mean of the)
responses at new feature values. Its input arguments should be out: output
produced by train.fun, and newx: feature values at which we want to make
predictions.}

\item{alpha}{Miscoverage level for the prediction intervals, i.e., intervals
with coverage 1-alpha are formed. Default for alpha is 0.1.}

\item{rho}{Fraction of observations to use in the training set. This is
rounded to the nearest integer below, and defaults to 0.5.}

\item{w}{Weights, in the case of covariate shift. This should be a vector of
length n+n0, giving the weights (i.e., ratio of test to training feature
densities), at each of n+n0 the training and test points. Default is NULL,
which means that we take all weights to be 1.}

\item{mad.train.fun}{A function to perform training on the absolute residuals
i.e., to produce an estimator of E(R|X) where R is the absolute residual
R = |Y - m(X)|, and m denotes the estimator produced by train.fun.
This is used to scale the conformal score, to produce a prediction interval
with varying local width. The input arguments to mad.train.fun should be
x: matrix of features, and y: vector of absolute residuals. The default for
mad.train.fun is NULL, which means that no training is done on the absolute
residuals, and the usual (unscaled) conformal score is used. Note that if
mad.train.fun is non-NULL, then so must be mad.predict.fun (see next).}

\item{mad.predict.fun}{A function to perform prediction for the (mean of the)
absolute residuals at new feature values. Its input arguments should be
out: output produced by mad.train.fun, and newx: feature values at which we
want to make predictions. The default for mad.predict.fun is NULL, which
means that no local scaling is done for the conformal score, i.e., the
usual (unscaled) conformal score is used.}

\item{split}{Indices that define the data-split to be used (i.e., the indices
define the first half of the data-split, on which the model is trained).
Default is NULL, in which case the split is chosen randomly.}

\item{seed}{Integer to be passed to set.seed before defining the random 
data-split to be used. Default is NULL, which effectively sets no seed.
If both split and seed are passed, the former takes priority and the latter
is ignored.}

\item{verbose}{Should intermediate progress be printed out? Default is FALSE.}
}
\value{
A list with the following components: pred, lo, up, fit, split. The
  first three are matrices of dimension n0 x m, and the fourth is a matrix of
  dimension n x m. Recall that n0 is the number of rows of x0, and m is the
  number of tuning parameter values internal to predict.fun. In a sense, each
  of the m columns really corresponds to a different prediction function;
  see details below. Hence, the rows of the matrices pred, lo, up give
  the predicted value, and lower and upper confidence limits (from split
  conformal inference), respectively, for the response at the n0 points
  given in x0. The rows of fit give the fitted values for the n points given
  in x. Finally, split contains the indices used for the first half of the
  data-split.
}
\description{
Compute prediction intervals using split conformal inference.
}
\details{
For concreteness, suppose that we want to use the predictions from
  forward stepwise regression at steps 1 through 5 in the path. In this case,
  there are m = 5 internal tuning parameter values to predict.fun, in the
  notation used above, and each of the returned matrices pred, lo, up, fit
  will have 5 columns (one for each step of the forward stepwise path).
  The code is structured in this way so that we may defined a single pair of
  functions train.fun and predict.fun, over a set of m = 5 tuning parameter
  values, instead of calling the conformal function separately m = 5 times.
}
\references{
See "Algorithmic Learning in a Random World" by Vovk, Gammerman,
  Shafer (2005) as the definitive reference for conformal prediction; see
  also "Distribution-Free Predictive Inference for Regression" by Lei,
  G'Sell, Rinaldo, Tibshirani, Wasserman (2018) for another description; and
  "Conformal Prediction Under Covariate Shift" by Barber, Candes, Ramdas,
  Tibshirani (2019) for the weighted extension.
}
\seealso{
\code{\link{conformal.pred}}, 
  \code{\link{conformal.pred.jack}}, \code{\link{conformal.pred.roo}}
}
