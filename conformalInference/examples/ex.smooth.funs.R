## Smoothing splines: some simple univariate examples

# Generate some example training data, clearly heteroskedastic
set.seed(33)
n = 1000
x = runif(n,0,2*pi)
y = sin(x) + x*pi/30*rnorm(n)

# Generate some example test data
n0 = 1000
x0 = runif(n,0,2*pi)
y0 = sin(x0) + x0*pi/30*rnorm(n0)

# Smoothing spline training and prediction functions, where the smoothing
# parameter is chosen via cross-validation
funs = smooth.spline.funs(cv=TRUE)

# Jackknife conformal for out-of-sample prediction intervals (valid on average
# over new x0 values), we can use the special leave-one-out fitting function
out.jack = conformal.pred.jack(x, y, x0, alpha=0.1,
  train.fun=funs$train, predict.fun=funs$predict, special.fun=funs$special)

plot(c(), c(), xlab="x0", ylab="y0", xlim=range(x0),
     ylim=range(c(y0,out.jack$lo,out.jack$up)), col="white",     
     main=paste0("Jackknife (out-of-sample) prediction intervals\n",
       sprintf("Average length: %0.3f",mean(out.jack$up-out.jack$lo))))
o = order(x0)
segments(x0[o], out.jack$lo[o], x0[o], out.jack$up[o], col="pink")
lines(x0[o], out.jack$pred[o], lwd=2, col="red")
points(x0, y0, col="gray50")

# Split conformal for out-of-sample prediction intervals (valid on average over
# new x0 values)
out.split = conformal.pred.split(x, y, x0, alpha=0.1,
  train.fun=funs$train, predict.fun=funs$predict)

plot(c(), c(), xlab="x0", ylab="y0", xlim=range(x0),
     ylim=range(c(y0,out.split$lo,out.split$up)), col="white",     
     main=paste0("Split conformal (out-of-sample) prediction intervals\n",
       sprintf("Average length: %0.3f",mean(out.split$up-out.split$lo))))
o = order(x0)
segments(x0[o], out.split$lo[o], x0[o], out.split$up[o], col="lightgreen")
lines(x0[o], out.split$pred[o], lwd=2, col="darkgreen")
points(x0, y0, col="gray50")

# Rank-one-out split conformal for in-sample prediction intervals (valid on 
# averaged over observed x values)
out.roo = conformal.pred.roo(x, y, alpha=0.1,
  train.fun=funs$train, predict.fun=funs$predict)

plot(c(), c(), xlab="x", ylab="y", xlim=range(x),
     ylim=range(c(y,out.roo$lo,out.roo$up)), col="white",
     main=paste0("ROO-conformal (in-sample) prediction intervals\n",
       sprintf("Average length: %0.3f",mean(out.roo$up-out.roo$lo))))
o = order(x)
segments(x[o], out.roo$lo[o], x[o], out.roo$up[o], col="lightblue")
lines(x[o], out.roo$fit.all[o], lwd=2, col="blue")
points(x, y, col="gray50")

## Adjust now for heteroskedasticity, by using mad.train.fun and mad.predict.fun

# Jackknife conformal is much slower to run, because we can't use the special
# leave-one-out fitting function anymore
out.jack.local = conformal.pred.jack(x, y, x0, alpha=0.1,
  train.fun=funs$train, predict.fun=funs$predict,
  mad.train.fun=funs$train, mad.predict.fun=funs$predict, verb=TRUE)

plot(c(), c(), xlab="x0", ylab="y0", xlim=range(x0),
     ylim=range(c(y0,out.jack.local$lo,out.jack.local$up)), col="white",     
     main=paste0("Jackknife (out-of-sample) prediction intervals\n",
       sprintf("Average length: %0.3f",
               mean(out.jack.local$up-out.jack.local$lo))))
o = order(x0)
segments(x0[o], out.jack.local$lo[o], x0[o], out.jack.local$up[o], col="pink")
lines(x0[o], out.jack.local$pred[o], lwd=2, col="red")
points(x0, y0, col="gray50")

# Split conformal, also using smooth splines (under CV) to train on residuals
out.split.local = conformal.pred.split(x, y, x0, alpha=0.1,
  train.fun=funs$train, predict.fun=funs$predict,
  mad.train.fun=funs$train, mad.predict.fun=funs$predict)

plot(c(), c(), xlab="x0", ylab="y0", xlim=range(x0), 
     ylim=range(c(y0,out.split.local$lo,out.split.local$up)), col="white",     
     main=paste0("Split conformal (out-of-sample) prediction intervals\n",
       sprintf("Average length: %0.3f",
               mean(out.split.local$up-out.split.local$lo))))
o = order(x0)
segments(x0[o], out.split.local$lo[o], x0[o], out.split.local$up[o],
         col="lightgreen")
lines(x0[o], out.split.local$pred[o], lwd=2, col="darkgreen")
points(x0, y0, col="gray50")

# Rank-one-out split conformal, again using smoothing splines (under CV) to
# train on residuals
out.roo.local = conformal.pred.roo(x, y, alpha=0.1,
  train.fun=funs$train, predict.fun=funs$predict,
  mad.train.fun=funs$train, mad.predict.fun=funs$predict)

plot(c(), c(), xlab="x", ylab="y", xlim=range(x),
     ylim=range(c(y,out.roo.local$lo,out.roo.local$up)), col="white",
     main=paste0("ROO-conformal (in-sample) prediction intervals\n",
       sprintf("Average length: %0.3f",
               mean(out.roo.local$up-out.roo.local$lo))))
o = order(x)
segments(x[o], out.roo.local$lo[o], x[o], out.roo.local$up[o], col="lightblue")
lines(x[o], out.roo.local$fit.all[o], lwd=2, col="blue")
points(x, y, col="gray50")

# The reason the local from rank-one-out method here don't look entirely smooth
# is because they're actually a mishmash of two separate split procedures on
# each half of the data set (so they'd look smooth over the points within each
# half, but not necessarily overall)
