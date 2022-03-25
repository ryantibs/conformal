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
  if(requireNamespace("future.apply", quietly=TRUE))
    return (future.apply::future_sapply(index,fun))
  else
    return (sapply(index,fun))
}


one.lapply= function(index, fun){
  if(requireNamespace("future.apply", quietly=TRUE))
    return (future.apply::future_lapply(index,fun))
  else
    return (lapply(index,fun))
}



#' Helper function to join each B prediction intervals in Multi Split
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

