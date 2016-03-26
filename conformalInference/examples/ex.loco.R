## Use lasso + CV to choose a model, then test variable importance in the
## selected model

# Generate some example training data
set.seed(33)
n = 200; p = 500; s = 5
x = matrix(rnorm(n*p),n,p)
beta = 2*c(rnorm(s),rep(0,p-s)) 
y = x %*% beta + rnorm(n)

# Lasso training and prediction functions, using cross-validation over 100
# values of lambda, with the 1se rule
funs = lasso.funs(nlambda=100,cv=TRUE,cv.rule="1se")

# Split-sample LOCO analysis
out.loco = loco(x, y, alpha=0.1, train.fun=funs$train, predict.fun=funs$predict,
  active.fun=funs$active.fun, verbose=TRUE)
out.loco

# Plot the Wilcoxon intervals
ylim = range(out.loco$inf.wilcox[[1]][,2:3])
J = length(out.loco$active[[1]])
plot(c(), c(), xlim=c(1,J), ylim=ylim, xaxt="n",
     main="LOCO analysis from lasso + CV model",
     xlab="Variable", ylab="Confidence interval")
axis(side=1, at=1:J, labels=FALSE)
for (j in 1:J) {
  axis(side=1, at=j, labels=out.loco$active[[1]][j], cex.axis=0.75,
       line=0.5*j%%2)
}
abline(h=0)
segments(1:J, out.loco$inf.wilcox[[1]][,2],
         1:J, out.loco$inf.wilcox[[1]][,3],
         col="red", lwd=2)
