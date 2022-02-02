# Better seq function
Seq = function(a, b, ...) {
  if (a <= b) return(seq(a,b,...))
  else return(integer(0))
}

# Truncate function
trunc = function(x, a, b) {
  return(ifelse(x>a,ifelse(x<b,2*x-a-b,b-a),a-b))
}

# Match duplicate row function
match.row = function(mat, a) {
  return(which(rowSums(abs(scale(mat,center=a,scale=F)))==0))
}

# Compute a coverage interval from a grid of values
grid.interval = function(pts, rvals, qvals, below=TRUE) {
  a = rvals <= qvals
  if (below) ii = which(a)
  else ii = which(!a)

  # If no rvals are on the correct side of qvals, then make
  # up a trivial interval
  if (length(ii)==0) return(list(lo=0,up=0))

  # Else take the lowest and highest points, form an interval
  i1 = min(ii)
  i2 = max(ii)
  return(list(lo=pts[i1],up=pts[i2]))
}

# Weighted quantile function
weighted.quantile = function(v, prob, w=NULL, sorted=FALSE) {
  if (is.null(w)) w = rep(1,length(v))
  if (!sorted) { o = order(v); v = v[o]; w = w[o] }
  i = which(cumsum(w/sum(w)) >= prob)
  if (length(i)==0) return(Inf) # Can happen with infinite weights
  else return(v[min(i)])
}

##Automatically switch between parallel or serial computation




one.sapply= function(index, fun){
  if("future.apply" %in% rownames(installed.packages()))
    return (future.apply::future_sapply(index,fun))
  else
    return (sapply(index,fun))
}


one.lapply= function(index, fun){
  if("future.apply" %in% rownames(installed.packages()))
    return (future.apply::future_lapply(index,fun))
  else
    return (lapply(index,fun))
}
