# Plot a few examples of LOCO analyses
library(conformalInference)

##### First the local type #####

# Generate some training data
set.seed(33)
n = 1000; p = 6
x = matrix(runif(n*p,-1,1),n,p)
f1 = function(u){sin((u+1)*pi)*(u<0)}
f2 = function(u){sin(u*pi)}
f3 = function(u){sin((u+1)*pi)*(u>0)}

mu = rowSums(cbind(f1(x[,1]),f2(x[,2]),f3(x[,3])))
y = mu + rnorm(n,0,0.1)

# Form training and prediction functions for a low-dimensional additive
# model, using the gam package
library(gam)
train.fun = function(x,y) {
  p = ncol(x)
  datafr = data.frame(x)
  str = paste("y ~",paste(paste0("s(X",1:p,",df=",5,")"),collapse="+"))
  return(gam(formula(str),data=datafr))
}
predict.fun = function(out,x0) {
  return(predict(out,data.frame(x0)))
}

# Use LOCO, with ROO conformal to form the prediction intervals for
# excess test error
out.loco = loco.roo(x, y, alpha=0.1, vars=0,
  train.fun=train.fun, predict.fun=predict.fun, verb=TRUE)

# Compute "c-values", the proportion of intervals that have a left endpoint
# less than or equal to 0
cvals = colMeans(out.loco$lo[,,1] <= 0)
cvals

mar = c(4.25,4.25,3.25,1)
w = 4
h = 4.5
xlab = "Location"
ylab = "Interval"

# For each variable in consideration, plot the intervals
for (j in 1:p) {
  pdf(file=sprintf("fig/loco.am.comp%i.pdf",j),w=w,h=h)
  par(mar=mar)
  plot(x[,j],x[,j],ylim=range(c(out.loco$lo[,j,1],out.loco$up[,j,1])),
       xlab=xlab,ylab=ylab,main=paste("Component",j),col=NA)
  cols = ifelse(out.loco$lo[,j,1] <= 0, 1, 3)
  segments(x[,j],out.loco$lo[,j,1],x[,j],out.loco$up[,j,1],col=cols)
  abline(h=0, lty=2, lwd=2, col=2)
  graphics.off()
  cat(sprintf("fig/loco.am.comp%i.pdf\n",j))
}
cat("\n")

##### Next the global type #####

# Generate some example training data
set.seed(33)
n = 200; p = 500; s = 5
x = matrix(rnorm(n*p),n,p)
beta = 2*c(rnorm(s),rep(0,p-s)) 
y = x %*% beta + rnorm(n)

# Lasso training and prediction functions, using cross-validation over 100
# values of lambda, with the 1se rule
my.lasso.funs = lasso.funs(nlambda=100,cv=TRUE,cv.rule="1se")

# Split-sample LOCO analysis
out.lasso = loco(x, y, alpha=0.1, train.fun=my.lasso.funs$train,
  predict.fun=my.lasso.funs$predict, active.fun=my.lasso.funs$active.fun,
  verbose=TRUE)

# Print out Wilcoxon results
print(out.lasso, test="wilcox")

# Plot the Wilcoxon intervals
mar = c(4.25,4.25,3.25,1)
w = 5
h = 5.5

pdf(file="fig/loco.lasso.pdf",w=w,h=h)
par(mar=mar)
ylim = range(out.lasso$inf.wilcox[[1]][,2:3])
J = length(out.lasso$active[[1]])
plot(c(), c(), xlim=c(1,J), ylim=ylim, xaxt="n",
     main="LOCO Analysis from Lasso + CV Model",
     xlab="Variable", ylab="Confidence interval")
axis(side=1, at=1:J, labels=FALSE)
for (j in 1:J) {
  axis(side=1, at=j, labels=out.lasso$active[[1]][j], cex.axis=0.75,
       line=0.5*j%%2)
}
abline(h=0)
segments(1:J, out.lasso$inf.wilcox[[1]][,2],
         1:J, out.lasso$inf.wilcox[[1]][,3],
         col="red", lwd=2)
graphics.off()
cat("\nfig/loco.lasso.pdf\n\n")

# Now repeat the analysis with forward stepwise, using cross-validation over
# 50 steps, with the 1se rule

# Forward stepwise training and prediction functions
my.step.funs = lars.funs(type="step",max.steps=50,cv=TRUE,cv.rule="1se")

# Split-sample LOCO analysis
out.step = loco(x, y, alpha=0.1, train.fun=my.step.funs$train,
  predict.fun=my.step.funs$predict, active.fun=my.step.funs$active.fun,
  verbose=TRUE)

# Print out Wilcoxon results
print(out.step, test="wilcox")

# Plot the Wilcoxon intervals
mar = c(4.25,4.25,3.25,1)
w = 5
h = 5.5

pdf(file="fig/loco.step.pdf",w=w,h=h)
par(mar=mar)
ylim = range(out.step$inf.wilcox[[1]][,2:3])
J = length(out.step$active[[1]])
plot(c(), c(), xlim=c(1,J), ylim=ylim, xaxt="n",
     main="LOCO Analysis from Stepwise + CV Model",
     xlab="Variable", ylab="Confidence interval")
axis(side=1, at=1:J, labels=FALSE)
for (j in 1:J) {
  axis(side=1, at=j, labels=out.step$active[[1]][j], cex.axis=0.75,
       line=0.5*j%%2)
}
abline(h=0)
segments(1:J, out.step$inf.wilcox[[1]][,2],
         1:J, out.step$inf.wilcox[[1]][,3],
         col="red", lwd=2)
graphics.off()
cat("\nfig/loco.step.pdf\n\n")

