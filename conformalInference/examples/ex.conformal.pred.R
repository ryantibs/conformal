## Lasso: use a fixed sequence of 100 lambda values

# Generate some example training data
set.seed(33)
n = 100; p = 120; s = 10
x = matrix(rnorm(n*p),n,p)
beta = c(rnorm(s),rep(0,p-s)) 
y = x %*% beta + rnorm(n)

# Generate some example test data
n0 = 50
x0 = matrix(rnorm(n0*p),n0,p)
y0 = x0 %*% beta + rnorm(n0)

# Grab a fixed lambda sequence from one call to glmnet, then
# define lasso training and prediction functions
if (!require("glmnet",quietly=TRUE)) {
  stop("Package glmnet not installed (required for this example)!")
}
out.gnet = glmnet(x,y,nlambda=100,lambda.min.ratio=1e-3)
lambda = out.gnet$lambda
funs = lasso.funs(lambda=lambda)

# Conformal inference, and jacknife and split conformal versions
out.conf = conformal.pred(x, y, x0, alpha=0.1,
  train.fun=funs$train, predict.fun=funs$predict, verb=TRUE)

out.jack = conformal.pred.jack(x, y, x0, alpha=0.1,
  train.fun=funs$train, predict.fun=funs$predict)

out.split = conformal.pred.split(x, y, x0, alpha=0.1,
  train.fun=funs$train, predict.fun=funs$predict)

y0.mat = matrix(rep(y0,ncol(out.conf$lo)),nrow=n0)
cov.conf = colMeans(out.conf$lo <= y0.mat & y0.mat <= out.conf$up)
len.conf = colMeans(out.conf$up - out.conf$lo)
err.conf = colMeans((y0.mat - out.conf$pred)^2)

cov.jack = colMeans(out.jack$lo <= y0.mat & y0.mat <= out.jack$up)
len.jack = colMeans(out.jack$up - out.jack$lo)
err.jack = colMeans((y0.mat - out.jack$pred)^2)

cov.split = colMeans(out.split$lo <= y0.mat & y0.mat <= out.split$up)
len.split = colMeans(out.split$up - out.split$lo)
err.split = colMeans((y0.mat - out.split$pred)^2)

# Compare to parametric intervals from oracle linear regression
lm.orac = lm(y~x[,1:s])
out.orac = predict(lm.orac,newdata=list(x=x0[,1:s]),
  interval="predict", level=0.9)

cov.orac = mean(out.orac[,"lwr"] <= y0 & y0 <= out.orac[,"upr"])
len.orac = mean(out.orac[,"upr"] - out.orac[,"lwr"])
err.orac = mean((y0 - out.orac[,"fit"])^2)
  
# Plot average coverage 
plot(log(lambda),cov.conf,type="l",ylim=c(0,1),
     xlab="log(lambda)",ylab="Avg coverage",
     main=paste0("Conformal + lasso (fixed lambda sequence):",
       "\nAverage coverage"))
points(log(lambda),cov.conf,pch=20)
lines(log(lambda),cov.jack,col=3)
points(log(lambda),cov.jack,pch=20,col=3)
lines(log(lambda),cov.split,col=4)
points(log(lambda),cov.split,pch=20,col=4)
abline(h=cov.orac,lty=2,col=2)
legend("bottomleft",col=c(1,3,4,2),lty=c(1,1,1,2),
       legend=c("Conformal","Jackknife conformal",
         "Split conformal","Oracle"))

# Plot average length
plot(log(lambda),len.conf,type="l",
     ylim=range(len.conf,len.orac),
     xlab="log(lambda)",ylab="Avg length",
     main=paste0("Conformal + lasso (fixed lambda sequence):",
       "\nAverage length"))
points(log(lambda),len.conf,pch=20)
lines(log(lambda),len.jack,col=3)
points(log(lambda),len.jack,pch=20,col=3)
lines(log(lambda),len.split,col=4)
points(log(lambda),len.split,pch=20,col=4)
abline(h=len.orac,lty=2,col=2)
legend("topleft",col=c(1,3,4,2),lty=c(1,1,1,2),
       legend=c("Conformal","Jackknife conformal",
         "Split conformal","Oracle"))

# Plot test error
plot(log(lambda),err.conf,type="l",
     ylim=range(err.conf,err.orac),
     xlab="log(lambda)",ylab="Test error",
     main=paste0("Conformal + lasso (fixed lambda sequence):",
       "\nTest error"))
points(log(lambda),err.conf,pch=20)
lines(log(lambda),err.jack,col=3)
points(log(lambda),err.jack,pch=20,col=3)
lines(log(lambda),err.split,col=4)
points(log(lambda),err.split,pch=20,col=4)
abline(h=err.orac,lty=2,col=2)
legend("topleft",col=c(1,3,4,2),lty=c(1,1,1,2),
       legend=c("Conformal","Jackknife conformal",
         "Split conformal","Oracle"))

