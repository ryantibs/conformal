# Plot an example comparison of usual to locally-weighted conformal, on
# heteroskedastic data 
library(conformalInference)

# Generate some example training data, clearly heteroskedastic
set.seed(33)
n = 1000
x = runif(n,0,2*pi)
y = sin(x) + x*pi/30*rnorm(n)

# Generate some example test data
n0 = 5000
x0 = seq(0,2*pi,length=n0)
y0 = sin(x0) + x0*pi/30*rnorm(n0)
o = order(x0)

# Smoothing spline training and prediction functions, where the smoothing
# parameter is chosen via cross-validation
funs = smooth.spline.funs(cv=TRUE)

# Split conformal, using smoothing splines
out.split = conformal.pred.split(x, y, x0, alpha=0.1, seed=10,
  train.fun=funs$train, predict.fun=funs$predict)
out.split.cov = out.split$lo <= y0 & y0 <= out.split$up
out.split.len = out.split$up-out.split$lo

pcol = "gray50"
mar = c(4.25,4.25,3.25,1)
w = 5
h = 5.5

pdf(file="fig/sin.split.pdf",w=w,h=h)
par(mar=mar)
plot(c(), c(), xlab="X", ylab="Y", xlim=range(x0),
     ylim=range(c(y0,out.split$lo,out.split$up)), col="white",     
     main=paste0("Split Conformal Prediction Intervals\n",
       sprintf("Coverage: %0.3f, Average Length: %0.3f",
               mean(out.split.cov),mean(out.split.len))))
polygon(c(x0[o],rev(x0[o])), c(out.split$lo[o],rev(out.split$up[o])),
        col="pink", border=NA)
lines(x0[o], out.split$pred[o], lwd=2, col="red")
points(x, y, col=pcol)
graphics.off()
cat("fig/sin.split.pdf\n")

# Split conformal, using smooth splines for both mean and residual
# estimation
out.split.cov = out.split$lo <= y0 & y0 <= out.split$up
out.split.len = out.split$up-out.split$lo
     
out.split.local = conformal.pred.split(x, y, x0, alpha=0.1, seed=10,
  train.fun=funs$train, predict.fun=funs$predict,
  mad.train.fun=funs$train, mad.predict.fun=funs$predict)
out.split.local.cov = out.split.local$lo <= y0 & y0 <= out.split.local$up
out.split.local.len = out.split.local$up-out.split.local$lo

pdf(file="fig/sin.split.local.pdf",w=w,h=h)
par(mar=mar)
plot(c(), c(), xlab="X", ylab="Y", xlim=range(x0),
     ylim=range(c(y0,out.split.local$lo,out.split.local$up)), col="white",     
     main=paste0("Locally-Weighted Split Conformal\n",
       sprintf("Coverage: %0.3f, Average Length: %0.3f",
               mean(out.split.local.cov),mean(out.split.local.len))))
polygon(c(x0[o],rev(x0[o])),
        c(out.split.local$lo[o],rev(out.split.local$up[o])),
        col="lightblue", border=NA)
lines(x0[o], out.split.local$pred[o], lwd=2, col="blue")
points(x, y, col=pcol)
graphics.off()
cat("fig/sin.split.local.pdf\n")

# Plot local coverage and length

# Smoothing spline
s.split.cov = smooth.spline(x0,as.numeric(out.split.cov),cv=TRUE)$y
s.split.local.cov = smooth.spline(x0,as.numeric(out.split.local.cov),cv=TRUE)$y
s.split.len = smooth.spline(x0,as.numeric(out.split.len),cv=TRUE)$y
s.split.local.len = smooth.spline(x0,as.numeric(out.split.local.len),cv=TRUE)$y

pdf(file="fig/sin.coverages.pdf",w=w,h=h)
par(mar=mar)
plot(c(), c(), xlab="X", ylab="Local coverage", xlim=range(x0),
     ylim=c(0,1), main="Local Coverage of Prediction Intervals")
lines(x0, s.split.cov, col="red")
lines(x0, s.split.local.cov, col="blue")
legend("bottomright", c("red","blue"), lty=c(1,1), col=c("red","blue"),
       legend=c("Usual","Locally-weighted"))
graphics.off()
cat("fig/sin.coverages.pdf\n")

pdf(file="fig/sin.lengths.pdf",w=w,h=h)
par(mar=mar)
plot(c(), c(), xlab="X", ylab="Local length", xlim=range(x0),
     ylim=range(c(s.split.len,s.split.local.len)),
     main="Local Length of Prediction Intervals")
lines(x0, s.split.len, col="red")
lines(x0, s.split.local.len, col="blue")
legend("bottomright", c("red","blue"), lty=c(1,1), col=c("red","blue"),
       legend=c("Usual","Locally-weighted"))
graphics.off()
cat("fig/sin.lengths.pdf\n")



