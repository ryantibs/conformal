#' Multi Split conformal prediction intervals.
#'
#' Compute prediction intervals using multisplit conformal inference.
#'
#' @param x Matrix of features, of dimension (say) n x p.
#' @param y Vector of responses, of length (say) n.
#' @param x0 Matrix of features, each row being a point at which we want to
#'   form a prediction interval, of dimension (say) n0 x p.
#' @param train.fun A function to perform model training, i.e., to produce an
#'   estimator of E(Y|X), the conditional expectation of the response variable
#'   Y given features X. Its input arguments should be x: matrix of features,
#'   and y: vector of responses.
#' @param predict.fun A function to perform prediction for the (mean of the)
#'   responses at new feature values. Its input arguments should be out: output
#'   produced by train.fun, and newx: feature values at which we want to make
#'   predictions.
#' @param alpha Miscoverage level for the prediction intervals, i.e., intervals
#'   with coverage 1-alpha are formed. Default for alpha is 0.1.
#' @param rho It is a vector of split proportions of length B. Default is NULL.
#' @param w Weights, in the case of covariate shift. This should be a vector of
#'   length n+n0, giving the weights (i.e., ratio of test to training feature
#'   densities), at each of n+n0 the training and test points. Default is NULL,
#'   which means that we take all weights to be 1.
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
#'   define the first half of the data-split, on which the model is trained).
#'   Default is NULL, in which case the split is chosen randomly.
#' @param seed Integer to be passed to set.seed before defining the random
#'   data-split to be used. Default is NULL, which effectively sets no seed.
#'   If both split and seed are passed, the former takes priority and the latter
#'   is ignored.
#' @param verbose Should intermediate progress be printed out? Default is FALSE.
#' @param B number of replications. Default is 50.
#' @param lambda Smoothing parameter. Default is 0.
#' @param tau It is a smoothing parameter, whose value affects the behavior of the
#'  function joining the B intervals:
#'     tau=1-1/B  Bonferroni intersection method
#'     tau=0 unadjusted intersection
#' Default is 1-(B+1)/(2*B).
#'
#' @return A list with the following components: lo, up. They are matrices of
#'  dimension n0 x m. Recall that n0 is the number of rows of x0, and m is the
#'   number of tuning parameter values internal to predict.fun. In a sense, each
#'   of the m columns really corresponds to a different prediction function;
#'   see details below.
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
#' @details This function is based on the package \code{\link{future.apply}} to
#'  perform parallelization. If this package is not installed, then the function
#'   will proceed serially.
#'
#' @seealso \code{\link{conformal.pred.split}}
#' @references "Multi Split Conformal Prediction" by Solari, Djordjilovic (2021).
#' @export conformal.pred.msplit



conformal.pred.msplit = function(x, y, x0,train.fun, predict.fun, alpha=0.1,
                                 rho = NULL, w=NULL, mad.train.fun=NULL, mad.predict.fun=NULL,
                                 split=NULL, seed=NULL,
                                 verbose=FALSE,
                                 B=50,
                                 lambda = 0,
                                 tau = 1-(B+1)/(2*B)) {


  x = as.matrix(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)
  x0 = matrix(x0,ncol=p)
  n0 = nrow(x0)
  
  ## Is the future.apply library installed?
  paral = "future.apply" %in% rownames(installed.packages())
  

  check.args(x=x,y=y,x0=x0,alpha=alpha,train.fun=train.fun,
             predict.fun=predict.fun,mad.train.fun=mad.train.fun,
             mad.predict.fun=mad.predict.fun)


  if (!is.null(seed)) set.seed(seed)
  else{
    seed=B
  }
  
  
  if (is.null(rho) || !is.vector(rho) || length(rho)!=B){
    print("rho value either not provided or with uncorrent format. 0.5 is set as
          drfault value in all the B replications!")
    rho = rep(0.5,B)
  }
  
  split.fit = is.matrix(split) && nrow(split) = B && ncol(split) > 0 && ncol(split) < n
  
  if(!split.fit && !is.null(split)){
    stop("The 'split' argument should either be NULL or a matrix of dimension B x n1,
         where n1 is the number of observation in the training set !  ")
  }
  

  loB <- upB <- matrix(NA,nrow=n0,ncol=B)
  tr <- tau*B + .001
  
  check.pos.num(lambda)
  check.num.01(tau)

  if(paral){
    #Define structures for parallelization
    future::plan(future::multisession)
    options(future.rng.onMisuse="ignore")
}

  ## Compute the B prediction intervals

    lo_up<- one.lapply(1:B, function(bbb) {

      out<-conformal.pred.split(x, y, x0, train.fun, predict.fun,
                                            alpha=alpha*(1-tau) + (alpha*lambda)/B, rho[bbb],
                                            w, mad.train.fun,
                                            mad.predict.fun,
                                            split[bbb,], seed=seed+bbb,
                                            verbose)


      return(cbind(t(out$lo),t(out$up)))
    })


  # Redefine the appropriate data structures

  Y_lo_up=do.call(rbind, lo_up)
  m=nrow(lo_up[[1]])
  index=rep(rep(1:m),B)
  list_mat=one.lapply(split(seq_along(index), index),
                  function(m, ind) m[ind,], m = Y_lo_up)[order(unique(index))]
  lo<-  up <- matrix(NA,nrow=n0,ncol=m)


  # Join the B prediction intervals

  for(ii in 1:m){

    loB=list_mat[[ii]][,1:n0]
    upB=list_mat[[ii]][,-(1:n0)]
    Y = rbind(loB,upB)

    finalInt = one.sapply(1:n0, function(kk) interval.build(Y[,kk],B,tr))

    lo[,ii]<-finalInt[1,]
    up[,ii]<-finalInt[2,]

  }
  
  if(paral){
    ## To avoid CRAN check errors
    ## R CMD check: make sure any open connections are closed afterward
    future::plan(future::sequential)

}

  return(list(lo=lo,up=up))
}




#' Helper function to join each B prediction intervals
#'
#' For each point in x0, it combines the simulated B intervals into one.
#'
#'
#' @param yyy column vector of B lower bounds and B upper bounds
#' @param B number of replications
#' @param tr truncation threshold for the algorithm



interval.build=function(yyy,B,tr){



  h=rep(1:0,each=B)

  o = order(yyy,2-h)

  ys <- yyy[o]
  hs <- h[o]

  count <- 0
  leftend <- 0
  lo<-up<-0


  for (j in 1:(2*B) ){
    if ( hs[j]==1 ) {
      count <- count + 1

      if ( count > tr && (count - 1) <= tr) {
        leftend <- ys[j]
      }

    }

    else {
      if ( count > tr && (count - 1) <= tr) {
        rightend <- ys[j]
        lo <- leftend
        up <- rightend
      }

      count <- count - 1
    }
  }




  return(c(lo,up))
}
