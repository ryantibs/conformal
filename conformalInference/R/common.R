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
grid.interval = function(pts, vals, threshold, right.side=TRUE,
  method=c("linear","floor","ceil")) {

  method = match.arg(method)
  n = length(pts)

  above = vals >= threshold
  if (right.side) ii = which(above)
  else ii = which(!above)

  # If no vals are on the correct side of the threshold, then
  # make up a trivial interval
  if (length(ii)==0) return(list(lo=0,up=0))

  i1 = min(ii)
  i2 = max(ii)
  
  if (method=="floor") {
    return(list(lo=pts[i1],up=pts[i2]))
  }
  else if (method=="ceil") {
    return(list(lo=pts[max(i1-1,1)],up=pts[min(i2+1,n)]))
  }
  else {
    if (i1==1) xl = pts[i1]
    else {
      slope = (vals[i1]-vals[i1-1])/(pts[i1]-pts[i1-1])
      xl = pts[i1-1] + (threshold-vals[i1-1])/slope
    }

    if (i2==n) xr = pts[i2]
    else {
      slope = (vals[i2+1]-vals[i2])/(pts[i2+1]-pts[i2])
      xr = pts[i2] + (threshold-vals[i2])/slope
    }
    
    return(list(lo=xl,up=xr))
  }
}

# Compute threshold for split conformal, at level alpha.
# Residuals must be sorted (and positive)
conformal.quantile = function(res, alpha) {
  n = length(res)
  if (ceiling((n+1)*alpha) <= 1) { return(Inf) }
  return(res[ceiling((n+1)*(1-alpha))])
}
