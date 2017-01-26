#' Features and responses generation.
#'
#' Generate a feature matrix x, and response vector y, following a specified
#' setup.
#'
#' @param n,p The dimensions, as in n x p, of the feature matrix x. Note that n
#'   is also the length of the response vector y.
#' @param x.dist One of "normal", "binom", "sn", "mix", indicating the
#'   distribution for the entries of x. In the first three cases, the columns
#'   of x are filled out with i.i.d. normal, standard binomial, or skewed normal
#'   (with skewness parameter alpha = 5) draws, respectively. The third option
#'   requiring the \code{\link{sn}} package.In the last case, "mix", each column
#'   of x is drawn with equal probability from the above three distributions.
#'   Default is "normal".
#' @param cor One of "none", "pair", "auto", "rand", indicating the type of
#'   correlation among columns of x. If "pair", then x is post-multiplied by
#'   the symmetric square root of the p x p matrix whose diagonals are 1 and
#'   off diagonals are rho (see parameter below). If "auto", then each column of
#'   x is taken to be a random convex combination of its old value and the
#'   previous k (see parameter) columns. If "rand", then each column of x is
#'   taken to be a random convex combination of its old value and a randomly
#'   k other columns. Default is "none".
#' @param rho If cor is "pair", then this specifies the pairwise correlation.
#'   Default is 0.5.
#' @param k If cor is "auto" or "rand", then this specifies the number of other
#'   columns that are used in the lagged or random convex combinations (see
#'   above). Default is 5.
#' @param standardize Should the columns of x be standardized (i.e., centered,
#'   and scaled to have variance 1)? Default is TRUE.
#' @param mean.fun One of "linear", "additive", or "rotated", indicating the
#'   functional dependence of the mean on the variables. If "additive", then the
#'   mean vector mu is constructed to be an additive function of B-splines of
#'   the variables, where each variable is expanded into df B-splines (with
#'   knots spaced evenly over the range of observed measurements of the
#'   variable). If "rotated", then mu is set to be randomly rotated version of
#'   such an additive function. Default is "linear".
#' @param m Number of B-splines per variable in the "additive" or "rotated" 
#'   cases for mean.fun; default is 4. Ignored when mean.fun is "linear".
#' @param sparsity Either "strong" or weak". In the former case, s (see below)
#'   randomly chosen variables of the p total contribute to mu, with
#'   coefficients randomly set to +1 or -1. In the latter case, still s variable
#'   contribute to mu, as before, but now the remaining p-s variables contribute
#'   in a weak sense, in that their coefficients are randomly set to +1 or -1
#'   times a successively increasing power of the weak.decay factor (see below).
#'   In the "additive" or "rotated" cases for mean.fun, the above should be
#'   interpreted in the appropriate groupwise sense, e.g., s groups of spline
#'   functions are randomly chosen, each group corresponding to the expansion of
#'   a single variable.
#' @param s Number indicating the level of sparsity. Note that the choice s = p
#'   would correspond to a fully dense model. Default is round(log(p)).
#' @param weak.decay Number controlling the decay of the weak coefficients, in
#'   the case sparsity is "weak", as their magnitudes get set to a successive
#'   power of weak.decay. Default is 0.5.
#' @param mu.seed Seed to be set before generation of the mean vector mu
#'   (ensures that the mean is the same, over multiple calls to this function).
#'   Default is NULL, which effectively sets no seed.
#' @param error.dist One of "normal", "unif", "t", "sn", indicating the error
#'   distribution, where the last option "sn" is for the skewed normal, and
#'   requires the \code{\link{sn}} package. Default is "normal".
#' @param sigma Standard deviation of errors. When the error distribution does
#'   not have a second moment (i.e., t-distribution with 1 or 2 df), then sigma
#'   represents the MAD. Default is 1.
#' @param bval Magnitude of nonzero coefficients in the model for mu as a linear
#'   function of x (or as a linear function of B-splines of x in the "additive"
#'   and "rotated" cases for mean.fun). Default is 3. This is overriden if the
#'   snr parameter is specified (see next).
#' @param snr The signal-to-noise ratio, specifically var(mu)/sigma^2, where
#'   sigma^2 is the variance of the errors (see previous). The mean will be
#'   scaled to meet the desired signal-to-noise ratio. 
#' @param sigma.type Either "const" or "var", indicating constant or varying
#'   error variance. In the latter case, the error standard deviation sigma is
#'   multiplied by 1 + abs(mu)/mean(abs(mu)) (this gives a vector of the same
#'   length as mu, i.e., a different standard deviation per component).
#' @param df If error.dist is "t", then this specifies the degrees of freedom of
#'   the t-distribution for the errors. Default is 3.
#' @param alpha If error.dist is "sn", then this specifies the skewness
#'   parameter (i.e., shape parameter) of the skewed normal distribution for the
#'   errors. The location and scale parameters are set so that the resulting
#'   distribution always has mean zero and variance 1. Note that alpha > 0 means
#'   right-skewed and alpha < 0 means left-skewed. Default is 5.
#' @param omitted.vars Should important (active) variables be omitted from x?
#'   Induces correlation between the x and the error terms.
#' 
#' @details In the case of an "additive" or "rotated" mean, this function relies
#'   on the packages \code{\link{plyr}} and \code{\link{splines}}. If these
#'   packages are not installed, then the function will abort.
#'
#' @return A list with the following components: x, y, and mu.
#'
#' @author Ryan Tibshirani
#' @references "Distribution-Free Predictive Inference for Regression" by 
#'   Jing Lei, Max G'Sell, Alessandro Rinaldo, Ryan Tibshirani, and Larry
#'   Wasserman, https://arxiv.org/pdf/1604.04173.pdf, 2016.
#' @examples ## See examples for sam.funs function
#' @export sim.xy

sim.xy = function(n, p, x.dist=c("normal","binom","sn","mix"),
  cor=c("none","pair","auto","rand"), rho=0.5, k=5, standardize=TRUE,
  mean.fun=c("linear","additive","rotated"), m=4, sparsity=c("strong","weak"),
  s=round(log(p)), weak.decay=0.5, mu.seed=NULL,
  error.dist=c("normal","unif","t","sn"), sigma=1, bval=3, snr=NULL,
  sigma.type=c("const","var"), df=3, alpha=5, omitted.vars=FALSE) {

  # Check arguments
  check.pos.int(n)
  check.pos.int(p)
  x.dist = match.arg(x.dist)
  cor = match.arg(cor)
  check.num.01(rho)
  check.pos.int(k)
  check.bool(standardize)
  mean.fun = match.arg(mean.fun)
  if (length(m) !=1 || !is.numeric(m) || m<2 || m>10 || round(m) != m) {
    stop("m must be an integer between 2 and 10")
  }
  sparsity = match.arg(sparsity)
  if (length(s) !=1 || !is.numeric(s) || s<1 || s>p || round(s) != s) {
    stop("s must be an integer between 1 and p")
  }
  check.num.01(weak.decay)
  error.dist = match.arg(error.dist)
  check.pos.num(sigma)
  check.pos.num(bval)
  if (!is.null(snr)) check.pos.num(snr)
  sigma.type = match.arg(sigma.type)
  check.pos.int(df)
  check.num(alpha)
  check.bool(omitted.vars)
  
  ####################
  # Generate x matrix
  if (x.dist=="normal") x = matrix(rnorm(n*p),n,p)
  if (x.dist=="binom") x = matrix(rbinom(n*p,1,0.5),n,p)
  if (x.dist=="sn") {
    if (!require(sn)) {
      stop("Package sn not installed (required here)!")
    }
    x = matrix(rsnorm(n*p,alpha=5),n,p)
  }
  if (x.dist=="mix") {
    if (!require(sn)) {
      stop("Package sn not installed (required here)!")
    }
    x1 = matrix(rnorm(n*p),n,p)
    x2 = matrix(rbinom(n*p,1,0.5),n,p)
    x3 = matrix(rsnorm(n*p,alpha=5),n,p)
    u = runif(p)
    i1 = u <= 1/3
    i2 = 1/3 < u & u <= 2/3
    i3 = u > 2/3
    x = matrix(0,n,p)
    x[,i1] = x1[,i1]
    x[,i2] = x2[,i2]
    x[,i3] = x3[,i3]
  }

  # Pairwise correlation
  if (cor=="pair") {
    b = (-2*sqrt(1-rho)+2*sqrt((1-rho)+p*rho))/(2*p)
    a = b + sqrt(1-rho)
    Sig.half = matrix(b,p,p)
    diag(Sig.half) = a
    x = x %*% Sig.half
  }

  # Auto-correlation
  else if (cor=="auto") {
    for (j in 1:p) {
      mat = x[,max(1,j-k):j,drop=F]
      wts = runif(ncol(mat),0,1)
      wts = wts/sum(wts)
      x[,j] = rowMeans(scale(mat,F,1/wts))
    }
  }

  # Rando-correlation
  else if (cor=="rand") {
    for (j in 1:p) {
      others = sample((1:p)[-j],k,replace=T)
      mat = x[,c(others,j),drop=F]
      wts = runif(ncol(mat),0,1)
      wts = wts/sum(wts)
      x[,j] = rowMeans(scale(mat,F,1/wts))
    }
  }

  # Standardize, if necessary
  if (standardize) x = scale(x,T,T)/sqrt(n-1)

  ####################
  # Generate error vector
  if (error.dist=="normal") e = rnorm(n)
  else if (error.dist=="unif") e = runif(n,-sqrt(12)/2,sqrt(12)/2)
  else if (error.dist=="t") e = rtfun(n,df)
  else if (error.dist=="sn") {
    if (!require(sn)) {
      stop("Package sn not installed (required here)!")
    }
    e = rsnorm(n,alpha)
  }
    
  ####################
  # Generate mu vector, set mu-specific seed
  if (!is.null(mu.seed)) set.seed(mu.seed)
  
  # Linear mean function
  if (mean.fun=="linear") {
    # Sample s variables uniformly at random, define true coefficients
    vars = sample(p,s)
    beta = rep(0,p)
    beta[vars] = sample(c(-1,1)*bval,s,replace=T)
    if (sparsity=="weak") {
      beta[-vars] = sample(c(-1,1)*bval,p-s,replace=T)*weak.decay^(1:(p-s)) 
    }
    mu = x %*% beta
  }
  
  # Additive or rotated mean function
  else {
    # Check for plyr
    if (!require(plyr)) {
      stop("Package plyr not installed (required here)!") 
    }
    # Check for splines
    if (!require(splines)) {
      stop("Package splines not installed (required here)!")
    }

    # Transform the columns of x into a B-spline expansion
    mat = do.call(cbind, alply(x, 2, function(x.col) bs(x.col,df=m) ))

    # N.B. Sometimes bs will have a problem when a particular column of x is
    # a (centered) binary variable, but not always. In that case, it'll return
    # NaNs ... fix this? Disallow spline transformations of binary variables?
    
    # Sample s variable groups uniformly at random, define true coefficients
    vars = sample(p,s)
    vars.mat = (rep(vars,each=m)-1)*m + 1:m
    beta = rep(0,p*m)
    beta[vars.mat] = sample(c(-1,1)*bval,s*m,replace=T)
    if (sparsity=="weak") {
      beta[-vars.mat] = sample(c(-1,1)*bval,(p-s)*m,replace=T) *
        weak.decay^(1:((p-s)*m))
    }
    mu = mat %*% beta
        
    # Rotate, if neccessary
    if (mean.fun=="rotated") mu = t(qr.Q(qr(matrix(rnorm(n^2),n,n)))) %*% mu
  }

  # Constant or varying error variance?
  if (sigma.type=="const") sx = rep(1,n)
  else sx = 1 + 2*abs(mu)^3/mean(abs(mu)^3)
  sigma.vec = sigma*sx
  
  # Set the signal-to-noise ratio, if necessary
  if (!is.null(snr)) mu = mu * sqrt(snr) * sqrt(mean(sigma.vec^2))/sd(mu)
  
  ####################
  # Finally, construct y
  y = mu + e*sigma.vec

  # In the omitted variables setup, we must delete half of the features; but 
  # make sure to delete half of the active ones, as well
  if (omitted.vars) {
    p1 = length(vars); p2 = p-p1
    x = cbind(x[,vars][,sample(p1,floor(p1/2))],
      x[,-vars][,sample(p2,ceiling(p2/2))])
  }
  
  return(list(x=x,y=y,mu=mu))
}

rtfun = function(n, df) {
  e = rt(n,df=df)
  if (df > 2) e = e / sqrt(df/(df-2))
  else e = e / mad(e)
  return(e)
}

rsnorm = function(n, alpha=5) {
  del = alpha/sqrt(1+alpha^2)
  xi = del*sqrt(2/pi)
  om = sqrt(1-2*del^2/pi)
  return((rsn(n,alpha=alpha)-xi)/om)
}

#' Master function for running simulations.
#'
#' Run a set of simulations with the specified configuration.
#'
#' @param n,p Number of observations and features, respectively.
#' @param conformal.pred.funs A list containing conformal inference functions
#'   to use. Each function here must take three arguments: x, y, x0 --- and no
#'   more. These will typically be defined with a wrapper around, e.g.,
#'   \code{\link{conformal.pred}} or \code{\link{conformal.pred.split}}. See
#'   examples below.
#' @param n0 The number of points at which to make predictions. Default is n.
#' @param in.sample Should the original x points be used for predictions (within
#'   each repetition)? If TRUE, then the n0 argument (above) is ignored. Default
#'   is FALSE.
#' @param nrep Number of repetitions over which to average results. Default is 
#'   20.
#' @param seed Seed to be set for the overall random number generation, i.e.,
#'   set before repetitions are begun (for reproducibility of the simulation
#'   results). Default is NULL, which effectively sets no seed.
#' @param verbose Should intermediate progress be printed out? Default is FALSE.
#' @param file,file.rep Name of a file to which simulation results will be saved
#'   (using saveRDS), and a number of repetitions after which intermediate
#'   results will be saved. Setting file to NULL is interpreted to mean that no
#'   simulations results should be saved; setting file.rep to 0 is interpreted
#'   to mean that simulations results should be saved at the very end, i.e., no
#'   intermediate saving. Defaults are NULL and 5, respectively.
#' @param x.dist,cor,rho,k,standardize,mean.fun,m,sparsity,s,weak.decay,
#'   error.dist,sigma,bval,snr,sigma.type,df,alpha,omitted.vars Arguments to 
#'   pass to \code{\link{sim.xy}}, see the latter's help file for details.
#'
#' @return A list with the following components.
#'   \itemize{
#'   \item ave.cov, ave.len, ave.err: lists of length m, where m denotes the
#'    number of conformal methods evaluated (i.e., the length of the list
#'    conformal.pred.funs). Each element of ave.cov corresponds to a conformal
#'    method considered, and is a vector whose length is equal to the number of
#'    tuning parameters internal to this method, containing the empirical
#'    coverages of conformal intervals, averaged across the repetitions. The
#'    lists ave.len and ave.err are analogous, except they contain the average
#'    lengths of conformal intervals, and average test errors, respectively.
#'   \item ave.opt, ave.tim: lists of length m, analogous to ave.cov, ave.len,
#'    and ave.err (see above), but where ave.opt contains the average (relative)
#'    optimism of the methods, a unitless measure of model complexity defined
#'    by (test error - training error) / (test error); and ave.tim contains
#'    the average runtimes (in seconds) of the methods.
#'   \item sd.cov, sd.len, sd.err, sd.opt, sd.tim: same as above, but containing
#'     deviations, rather than averages, across the reptitions.
#'   \item cov, len, err, opt, tim: lists of length m, reporting the full set
#'     of metrics across repetitions (rather than averages or standard
#'     deviations).
#'   \item ave.best.err, ave.best.len: matrices of dimension 5 x m, providing a
#'    summary of the average coverage, length, error, optimism, and time
#'    metrics (see description above), but at only one particular tuning
#'    parameter value for each method considered --- this is the tuning parameter
#'    value that either gives the best average error, or best average length.
#'   \item sd.best.err, sd.best.len: matrices of dimension 5 x m,
#'    providing a summary of the standard deviation of coverage, length, error,
#'    optimism, and time metrics for each method's best tuning parameter value
#'    (analogous to the construction of ave.best.err and ave.best.len).
#'   }
#'
#' @seealso \code{\link{sim.xy}}
#' @author Ryan Tibshirani
#' @references "Distribution-Free Predictive Inference for Regression" by 
#'   Jing Lei, Max G'Sell, Alessandro Rinaldo, Ryan Tibshirani, and Larry
#'   Wasserman, https://arxiv.org/pdf/1604.04173.pdf, 2016.
#' @example examples/ex.sim.master.R
#' @export sim.master

sim.master = function(n, p, conformal.pred.funs, n0=n, in.sample=FALSE, nrep=20,
  seed=NULL, verbose=FALSE, file=NULL, file.rep=5, 
  x.dist=c("normal","binom","sn","mix"), cor=c("none","pair","auto","rand"),
  rho=0.5, k=5, standardize=TRUE, mean.fun=c("linear","additive","rotated"),
  m=4, sparsity=c("strong","weak"), s=round(log(p)), weak.decay=0.5,
  mu.seed=NULL, error.dist=c("normal","unif","t","sn"), sigma=1, bval=3,
  snr=NULL, sigma.type=c("const","var"), df=3, alpha=5, omitted.vars=FALSE) {
  
  this.call = match.call()
  
  # Check arguments
  check.pos.int(n)
  check.pos.int(p)
  if (is.null(conformal.pred.funs)) {
    stop(paste("Must provide at least one method for",
               "conformal inference, through conformal.pred.funs"))
  }
  if (!is.list(conformal.pred.funs)
      || !all(sapply(conformal.pred.funs,is.function))) {
    stop("conformal.pred.funs must be a list of functions")
  }
  check.pos.int(n0)
  check.bool(in.sample)
  check.pos.int(nrep)
  check.int(file.rep)
  check.bool(omitted.vars)
  if (!is.null(seed)) set.seed(seed)
  
  N = length(conformal.pred.funs)
  nms = names(conformal.pred.funs)
  if (is.null(nms)) nms = paste("Method",1:N)
  
  cov = len = err = opt = tim = vector(mode="list",length=N)
  names(cov) = names(len) = names(err) = names(opt) = names(tim) = nms
  for (j in 1:N) {
    cov[[j]] = len[[j]] = err[[j]] = opt[[j]] = matrix(NA,nrep,1)
    tim[[j]] = rep(NA,nrep)
  }
  filled = rep(FALSE,N)

  # Loop through the repetitions
  for (i in 1:nrep) {
    if (verbose) {
      cat(sprintf("Simulation %i (of %i) ...\n",i,nrep))
      cat("\tGenerating data ...\n")
    }
    
    # Generate x, y, x0, y0 --- all together
    if (in.sample) { x0 = x; n0 = n }
    obj = sim.xy(n+n0,p,x.dist=x.dist,cor=cor,rho=rho,k=k,
      standardize=standardize,mean.fun=mean.fun,m=m,sparsity=sparsity,s=s,
      weak.decay=weak.decay,mu.seed=mu.seed,error.dist=error.dist,sigma=sigma,
      bval=bval,snr=snr,sigma.type=sigma.type,df=df,omitted.vars=omitted.vars)
    x = obj$x[1:n,]
    y = obj$y[1:n]
    x0 = obj$x[(n+1):(n+n0),]
    y0 = obj$y[(n+1):(n+n0)]
    
    # Loop through the conformal methods
    for (j in 1:N) {
      if (verbose) {
        cat(sprintf("\tApplying conformal inference method %i (of %i) ...\n",
                    j,N))
      }
      
      tt = system.time({
        tryCatch({
          # Apply the conformal method in hand
          out = conformal.pred.funs[[j]](x=x,y=y,x0=x0)
          y.mat = matrix(rep(y,ncol(out$lo)),nrow=n)
          y0.mat = matrix(rep(y0,ncol(out$lo)),nrow=n0)
   
          # Populate empty matrices for our metrics, of appropriate dimension
          if (!filled[j]) {
            cov[[j]] = len[[j]] = err[[j]] = opt[[j]] =
              matrix(NA,nrep,ncol(out$lo))
            filled[j] = TRUE
            # N.B. Filling with NAs is important, because the filled flag could
            # be false for two reasons: i) we are at the first iteration, or ii)
            # we've failed in all previous iters to run the conformal method
          }

          # Record all of our metrics
          cov[[j]][i,] = colMeans(out$lo <= y0.mat & y0.mat <= out$up)
          len[[j]][i,] = colMeans(out$up - out$lo)
          train.err = colMeans((out$fit - y.mat)^2)
          test.err = colMeans((out$pred - y0.mat)^2)
          err[[j]][i,] = test.err
          opt[[j]][i,] = (test.err - train.err) / test.err
        }, error = function(err) {
          if (verbose) {
            cat(paste("\t\tOops! Something went wrong, see error message",
                      "below; recording all metrics here as NAs ...\n"))
            cat("\t\t***** Error message *****\n")
            cat(sprintf("\t\t%s\n",err$message))
            cat("\t\t*** End error message ***\n")
          }
          # N.B. No need to do anything, the metrics are already filled with NAs
        })
      })[1]

      tim[[j]][i] = tt
    }

    # Save intermediate results?
    if (!is.null(file) && file.rep > 0 && i %% file.rep == 0) {
      saveRDS(list(cov,len,err,opt,tim), file)
    }
  }

  # Average out the various metrics over the simulations
  ave.cov = ave.len = ave.err = ave.opt = ave.tim = vector(mode="list",length=N)  
  sd.cov = sd.len = sd.err = sd.opt = sd.tim = vector(mode="list",length=N)
  names(ave.cov) = names(ave.len) = names(ave.err) = nms 
  names(sd.cov) = names(sd.len) = names(sd.err) = nms
  names(ave.opt) = names(ave.tim) = names(sd.opt) = names(sd.tim) = nms
  for (j in 1:N) {
    ave.cov[[j]] = colMeans(cov[[j]],na.rm=T)
    ave.len[[j]] = colMeans(len[[j]],na.rm=T)
    ave.err[[j]] = colMeans(err[[j]],na.rm=T)
    ave.opt[[j]] = colMeans(opt[[j]],na.rm=T)
    ave.tim[[j]] = mean(tim[[j]],na.rm=T)
    sd.cov[[j]] = apply(cov[[j]],2,sd,na.rm=T)
    sd.len[[j]] = apply(len[[j]],2,sd,na.rm=T)
    sd.err[[j]] = apply(err[[j]],2,sd,na.rm=T)
    sd.opt[[j]] = apply(opt[[j]],2,sd,na.rm=T)
    sd.tim[[j]] = sd(tim[[j]],na.rm=T)
  }

  # Populate summary tables (highlighting best error or length for each method) 
  ave.best.err = ave.best.len =  sd.best.err = sd.best.len = matrix(0,5,N)
  for (j in 1:N) {
    i.err = which.min(ave.err[[j]])
    i.len = which.min(ave.len[[j]])
    ave.best.err[,j] = c(ave.cov[[j]][i.err],ave.len[[j]][i.err],
                  ave.err[[j]][i.err],ave.opt[[j]][i.err],ave.tim[[j]])
    ave.best.len[,j] = c(ave.cov[[j]][i.len],ave.len[[j]][i.len],
                  ave.err[[j]][i.len],ave.opt[[j]][i.len],ave.tim[[j]])
    sd.best.err[,j] = c(sd.cov[[j]][i.err],sd.len[[j]][i.err],
                 sd.err[[j]][i.err],sd.opt[[j]][i.err],sd.tim[[j]])
    sd.best.len[,j] = c(sd.cov[[j]][i.len],sd.len[[j]][i.len],
                 sd.err[[j]][i.len],sd.opt[[j]][i.len],sd.tim[[j]])
  }
  colnames(ave.best.err) = colnames(ave.best.len) = nms
  rownames(ave.best.err) = rownames(ave.best.len) =
    c("Coverage","Length","Test error","Optimism","Time")
  colnames(sd.best.err) = colnames(sd.best.len) = nms
  rownames(sd.best.err) = rownames(sd.best.len) =
    c("Coverage","Length","Test error","Optimism","Time")
    
  out = list(ave.cov=ave.cov,ave.len=ave.len,ave.err=ave.err,
    ave.opt=ave.opt,ave.tim=ave.tim,
    sd.cov=sd.cov,sd.len=sd.len,sd.err=sd.err,
    sd.opt=sd.opt,sd.tim=sd.tim,
    cov=cov,len=len,err=err,opt=opt,tim=tim,
    ave.best.err=ave.best.err,ave.best.len=ave.best.len,
    sd.best.err=sd.best.err,sd.best.len=sd.best.len,
    call=this.call)
  class(out) = "sim"
  
  # Save results?
  if (!is.null(file)) saveRDS(out, file)

  return(out)
}

#' Print function for sim object.
#'
#' Summarize and print the results of a set of simulations, stored an object
#' of class sim (produced by sim.master).
#'
#' @param x The sim object.
#' @param type One of "all", "err", or "len" describing the type of summary to 
#'   be used. If "err", then for each method, the tuning parameter value that
#'   gives the best average test error is selected, and the coverage and length
#'   properties at this tuning parameter value are report. If "len", then the
#'   analogous thing is done, but for the tuning parameter value that gives the
#'   best average length. If "all", then the average coverage, length, and test
#'   error results are shown across all tuning parameter values. Default is
#'   "err".
#' @param se Should standard errors be displayed (in parantheses)? For each
#'   metric, this is defined as the standard deviation divided by the square
#'   root of the number of simulation repetitions. Default is TRUE.
#' @param digits Number of digits to display. Default is 3. 
#' @param ... Other arguments (currently not used).
#'
#' @export 

print.sim = function(x, type=c("err","len","all"), se=TRUE, digits=3, ...) {
  type = match.arg(type)
  
  if (!is.null(x$call)) {
    cat("\nCall:\n")
    dput(x$call)
  }
  nrep = nrow(x$err[[1]])
  
  if (type=="err") {
    cat(paste0("\nSummary (displaying results, for each method, at the tuning",
               "\nparameter value that minimizes test error):\n\n"))
    tab = round(x$ave.best.err,digits)
    if (se) {
      tab[1:3,] = matrix(paste0(tab[1:3,]," (",
           round(x$sd.best.err[1:3,]/sqrt(nrep),digits),")"),
           ncol=ncol(tab))
    }
    print(tab,quote=F)
  }
  else if (type=="len") {
    cat(paste0("\nSummary (displaying results, for each method, at the tuning",
              "\nparameter value that minimizes length):\n\n"))
    tab = round(x$ave.best.len,digits)
    if (se) {
      tab[1:3,] = matrix(paste0(tab[1:3,]," (",
           round(x$sd.best.len[1:3,]/sqrt(nrep),digits),")"),
           ncol=ncol(tab))
    }
    print(tab,quote=F)
  }
  else {
    cat("\nSummary (over all tuning parameter values):\n")
    N = length(x$ave.cov)
    nms = names(x$ave.cov)
    if (is.null(nms)) nms = paste("Method",1:N)
    
    for (j in 1:N) {
      tab = round(rbind(x$ave.cov[[j]],
        x$ave.len[[j]],x$ave.err[[j]]),digits=3)
      
      if (se) {
        tab = matrix(paste0(tab," (",
          round(rbind(x$sd.cov[[j]],x$sd.len[[j]],
                      x$sd.err[[j]]),digits=3),
          ")"),nrow(tab),ncol(tab))
      }
      
      rownames(tab) = c("Coverage","Length","Test error")
      colnames(tab) = paste0(c("Opt ",rep("",length(x$ave.opt[[j]])-1)),
                round(x$ave.opt[[j]],digits))

      cat(paste("\n-------",nms[j],"-------\n"))
      cat(paste("Average time:",round(x$ave.tim[[j]],digits),"seconds\n"))
      print(tab,quote=F)
    }
  }

  cat("\n")
  invisible()
}

#' Plot function for sim object.
#'
#' Plot the results of a set of simulations, stored in an object of class sim
#' (produced by sim.master).
#'
#' @export

plot.sim = function(x, method.nums=1:length(x$ave.cov), method.names=NULL, 
  cols=1:8, log="", xlim=NULL, se=FALSE, main=NULL, cex.main=1.25,
  xlab="Relative optimism", legend.pos=c("bottomleft","",""),
  make.pdf=FALSE, fig.dir=".", file.prefix="", h=6, w=6, mar=NULL, ...) {
 
  if (is.null(method.names)) method.names = paste("Method",method.nums)
  if (is.null(mar) && is.null(main)) mar = c(4.25,4.25,1,1)
  if (is.null(mar) && !is.null(main)) mar = c(4.25,4.25,2.25,1)
  if (is.null(main)) main = ""
  cols = rep(cols,length=length(method.nums))
  main = rep(main,length=3)
  legend.pos = rep(legend.pos,length=3)
  ii = method.nums
  nrep = nrow(x$cov[[1]])

  # Coverage plot
  plot.one(x$ave.opt[ii], x$ave.cov[ii],
           se.y=lapply(x$sd.cov[ii], function(v) return(v/sqrt(nrep))), se=se,
           names=method.names, cols=cols, log=log, xlim=xlim, ylim=c(0,1),
           main=main[1], legend.pos=legend.pos[1], make.pdf=make.pdf,
           fig.dir=fig.dir, file.prefix=paste0(file.prefix,".cov"), h=h, w=w,
           mar=mar, xlab=xlab, ylab="Coverage", trunc=F, ...)

  # Length plot
  plot.one(x$ave.opt[ii], x$ave.len[ii], 
           se.y=lapply(x$sd.len[ii], function(v) return(v/sqrt(nrep))), se=se,
           names=method.names, cols=cols, log=log, xlim=xlim, ylim=NULL,
           main=main[2], legend.pos=legend.pos[2], make.pdf=make.pdf,
           fig.dir=fig.dir, file.prefix=paste0(file.prefix,".len"), h=h, w=w,
           mar=mar, xlab=xlab, ylab="Length", trunc=T, ...)

  # Test error plot
  plot.one(x=x$ave.opt[ii], y=x$ave.err[ii],
           se.y=lapply(x$sd.err[ii], function(v) return(v/sqrt(nrep))), se=se,
           names=method.names, cols=cols, log=log, xlim=xlim, ylim=NULL,
           main=main[3], legend.pos=legend.pos[3], make.pdf=make.pdf,
           fig.dir=fig.dir, file.prefix=paste0(file.prefix,".err"), h=h, w=w,
           mar=mar, xlab=xlab, ylab="Test error", trunc=T, ...)

}

plot.one = function(x, y, se.y, se=FALSE, names=NULL, cols=1:8, log="",
  xlim=NULL, ylim=NULL, main=NULL, cex.main=1.25, legend.pos="bottomleft",
  make.pdf=FALSE, fig.dir=".", file.prefix="", h=6, w=6, mar=NULL,
  xlab=NULL, ylab=NULL, trunc=FALSE, ...) {

  N = length(x); cols = rep(cols,length=N)
  xall = unlist(x)
  if (is.null(xlim)) xlim = range(xall,finite=T)
  trans = ifelse(xlim[1] <= 0, -xlim[1] + 1, 0)
  xseq = seq(xlim[1],xlim[2],length=10)
  boo = ind = vector(mode="list",length=N)
  for (j in 1:N) {
    boo[[j]] = xlim[1] <= x[[j]] & x[[j]] <= xlim[2]
    ind[[j]] = which(boo[[j]])
  }
  
  if (make.pdf) pdf(file=sprintf("%s/%s.pdf",fig.dir,file.prefix),w,h)
  else dev.new()
  par(mar=mar)

  if (is.null(ylim)) {
    if (se) ylim = range(c((unlist(y)-unlist(se.y))[unlist(boo)],
              (unlist(y)+unlist(se.y))[unlist(boo)]),finite=T)
    else ylim = range(unlist(y)[unlist(boo)],finite=T)
    if (trunc) ylim = pmax(ylim,0)
  }
  
  xpts = x[[1]][ind[[1]]]+trans
  plot(xpts, y[[1]][ind[[1]]], xaxt=ifelse(trans==0,"s","n"),
       xlim=xlim+trans, ylim=ylim, type="l", col=cols[1], log=log,
       xlab=xlab, ylab=ylab, main=main, cex.main=cex.main)
  if (trans!=0) axis(side=1, at=xseq+trans, labels=round(xseq,1))
  points(xpts, y[[1]][ind[[1]]], pch=20, col=cols[1])
  if (se) segments(xpts, (y[[1]]-se.y[[1]])[ind[[1]]],
                   xpts, (y[[1]]+se.y[[1]])[ind[[1]]],
                   lwd=0.5, col=cols[1])        
  for (j in Seq(2,N)) {
    xpts = x[[j]][ind[[j]]]+trans
    lines(xpts, y[[j]][ind[[j]]], col=cols[j])
    points(xpts, y[[j]][ind[[j]]], pch=20, col=cols[j])
    if (se) segments(xpts, (y[[j]]-se.y[[j]])[ind[[j]]],
                     xpts, (y[[j]]+se.y[[j]])[ind[[j]]],
                     lwd=0.5, col=cols[j])   
  }
  if (legend.pos != "") legend(legend.pos,legend=names, col=cols,lty=1)
  if (make.pdf) graphics.off()
}

#' Print function for latex-style tables.
#'
#' Print a given table in format digestable by latex.
#' 
#' @export print.tex

print.tex = function(tab, tab.se=NULL, digits=3, file=NULL, align="l") {
  tab = round(tab,digits)
  n = nrow(tab); m = ncol(tab)
  rownms = rownames(tab)
  colnms = colnames(tab)
  if (!is.null(tab.se)) {
    tab = matrix(paste0(tab, " (", round(tab.se,digits), ")"), ncol=m)
  }
  
  if (is.null(file)) file = ""
  cat("", file=file, append=F) # Wipe the file, or the console
  cat(paste0("\\begin{tabular}{|",paste0(rep(align,m+1),collapse="|"),"|}\n"),
      file=file, append=T)
  cat("\\hline\n", file=file, append=T)
  if (!is.null(colnms)) {
    for (j in 1:m) cat(paste0("& ", colnms[j]," "), file=file, append=T)
    cat("\\\\\n", file=file, append=T)
    cat("\\hline\n", file=file, append=T)
  }
  for (i in 1:n) {
    cat(paste0(rownms[i], " "), file=file, append=T)
    for (j in 1:m) cat(paste0("& ",tab[i,j]," "), file=file, append=T)
    cat("\\\\\n", file=file, append=T)
    cat("\\hline\n", file=file, append=T)
  }
  cat(paste0("\\end{tabular}\n"), file=file, append=T)
}
