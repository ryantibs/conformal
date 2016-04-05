## Use lasso + CV to choose a model, then test variable importance in the
## selected model

# Generate some example training data
set.seed(11)
n = 200; p = 500; s = 10
x = matrix(rnorm(n*p),n,p)
beta = 8*c(rnorm(s),rep(0,p-s)) 
y = x %*% beta + rnorm(n)

# Lasso training and prediction functions, using cross-validation over 100
# values of lambda, with the 1se rule
funs = lasso.funs(nlambda=100,cv=TRUE,cv.rule="1se")

# In-sample rank-one-out conformal inference 
out.roo = conformal.pred.roo(x, y, alpha=0.1,
  train.fun=funs$train, predict.fun=funs$predict)

# Look at which variables were nonzero in the lasso + CV model fit to the whole
# data set
out.all = out.roo$out.all
vars = which(coef(out.all$glmnet.fit, s=out.all$lambda.1se) != 0)

# Look at error inflation due to variable dropping, among vars
out.loco = loco.roo(x, y, alpha=0.1, vars=vars,
  train.fun=funs$train, predict.fun=funs$predict, verbose=TRUE)

# Compute "c-values", the proportion of intervals that have a left endpoint
# less than or equal to 0
cvals = colMeans(out.loco$lo[,,1] <= 0)
names(cvals) = vars
cvals

xlab = "Sample number"
ylab = "Interval"

# For each dropped variable in consideration, plot the intervals
for (j in 1:length(vars)) {
  plot(c(),c(),xlim=c(1,n),ylim=range(c(out.loco$lo[,j,1],out.loco$up[,j,1])),
       xlab=xlab,ylab=ylab,main=sprintf("Variable %i, c-value = %0.3f",
                             vars[j], cvals[j]))
  cols = ifelse(out.loco$lo[,j,1] <= 0, 1, 3)
  segments(1:n,out.loco$lo[,j,1],1:n,out.loco$up[,j,1],col=cols)
  abline(h=0, lty=2, lwd=2, col=2)
}
