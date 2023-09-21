check.args = function(x=NULL, y=NULL, x0=NULL, alpha=NULL,
  train.fun=NULL, predict.fun=NULL, mad.train.fun=NULL, mad.predict.fun=NULL,
  special.fun=NULL) {
  
  if (is.null(x) || !is.numeric(x)) stop("x must be a numeric matrix")
  if (is.null(y) || !is.numeric(y)) stop("y must be a numeric vector")
  if (nrow(x) != length(y)) stop("nrow(x) and length(y) must match")
  if (is.null(x0) || !is.numeric(x0)) stop("x0 must be a numeric matrix")
  if (ncol(x) != ncol(x0)) stop("ncol(x) and ncol(x0) must match")
  check.num.01(alpha)
  if (is.null(train.fun) || !is.function(train.fun))
    stop("train.fun must be a function")
  if (is.null(predict.fun) || !is.function(predict.fun))
    stop("predict.fun must be a function")
  if (!is.null(mad.train.fun) && !is.function(mad.train.fun)) 
    stop("mad.train.fun must be a function")
  if (!is.null(mad.predict.fun) && !is.function(mad.predict.fun)) 
    stop("mad.predict.fun must be a function")
  if ((!is.null(mad.train.fun) && is.null(mad.predict.fun)) ||
      (is.null(mad.train.fun) && !is.null(mad.predict.fun)))
    stop("mad.train.fun and mad.predict.fun must both be provided")
  if (!is.null(special.fun) && !is.function(special.fun)) 
    stop("special.fun must be a function")
}

check.args.surv = function(x=NULL, t=NULL, d=NULL, tau=NULL, x0=NULL, w=NULL,
                           cens.model=NULL, alpha=NULL, train.fun=NULL, predict.fun=NULL, 
                           mad.train.fun=NULL, mad.predict.fun=NULL, special.fun=NULL) {
  
  if (is.null(x) || !is.numeric(x)) stop("x must be a numeric matrix")
  if (is.null(t) || !is.numeric(t)) stop("t must be a numeric vector")
  if (nrow(x) != length(t)) stop("nrow(x) and length(t) must match")
  if (is.null(d) || !is.numeric(d)) stop("d must be a binary (0/1) vector")
  if (length(d[d!=0 & d!=1]) != 0) stop("d must be a binary (0/1) vector")
  if (length(d) != length(t)) stop("length(d) and length(t) must match")
  if (!is.null(w) && !is.numeric(w)) stop("w must be NULL or a numeric vector")
  if (!is.null(w) && length(w) != length(t)) 
    stop("length(w) and length(t) must match")
  if (!(cens.model %in% c("km","rsf","cox"))) 
    stop("cens.model must be km, rsf or cox")
  check.pos.num(tau)
  check.num.01(alpha)
  check.args = function(x, t, x0, alpha, train.fun,
    predict.fun, mad.train.fun, mad.predict.fun, special.fun)
}

check.bool = function(b) {
  if (is.null(b) || length(b)!=1 || !is.logical(b))
    stop(paste(deparse(substitute(b)),"must be a Boolean"))
}

check.num = function(a) {
  if (is.null(a) || length(a)!= 1 || !is.numeric(a))
    stop(paste(deparse(substitute(a)),"must be a number"))
}

check.int = function(i) {
  if (is.null(i) || length(i)!= 1 || !is.numeric(i) || round(i) != i)
    stop(paste(deparse(substitute(i)),"must be an integer"))
}

check.pos.num = function(a) {
  if (is.null(a) || length(a)!= 1 || !is.numeric(a) || a<0)
    stop(paste(deparse(substitute(a)),"must be a positive number"))
}

check.pos.int = function(i) {
  if (is.null(i) || length(i)!= 1 || !is.numeric(i) || round(i) != i || i<1)
    stop(paste(deparse(substitute(i)),"must be a positive integer"))
}

check.num.01 = function(a) {
  if (is.null(a) || length(a)!= 1 || !is.numeric(a) || a<0 || a>1)
    stop(paste(deparse(substitute(a)),"must be a number between 0 and 1"))
}