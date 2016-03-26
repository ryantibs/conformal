# Plot a few example of loco analyses
library(conformalInference)

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
