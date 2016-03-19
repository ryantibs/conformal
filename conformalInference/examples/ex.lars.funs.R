## Forward stepwise: use 30 path steps

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

# Stepwise training and prediction functions
funs = lars.funs(type="stepwise",max.steps=30)

# Split conformal inference 
out.split = conformal.pred.split(x, y, x0, alpha=0.1,
  train.fun=funs$train, predict.fun=funs$predict)

y0.mat = matrix(rep(y0,ncol(out.split$lo)),nrow=n0)
cov.split = colMeans(out.split$lo <= y0.mat & y0.mat <= out.split$up)
len.split = colMeans(out.split$up - out.split$lo)
err.split = colMeans((y0.mat - out.split$pred)^2)

# Compare to parametric intervals from oracle linear regression
lm.orac = lm(y~x[,1:s])
out.orac = predict(lm.orac, list(x=x0[,1:s]),
  interval="predict", level=0.9)

cov.orac = mean(out.orac[,"lwr"] <= y0 & y0 <= out.orac[,"upr"])
len.orac = mean(out.orac[,"upr"] - out.orac[,"lwr"])
err.orac = mean((y0 - out.orac[,"fit"])^2)
  
# Plot average coverage 
plot(1:length(cov.split),cov.split,type="l",ylim=c(0,1),
     xlab="Number of path steps",ylab="Avg coverage",
     main=paste0("Split conformal + stepwise:\nAverage coverage"))
points(1:length(cov.split),cov.split,pch=20)
abline(h=cov.orac,lty=2,col=2)
legend("bottomright",col=1:2,lty=1:2,
       legend=c("Split conformal","Oracle"))

# Plot average length
plot(1:length(len.split),len.split,type="l",
     ylim=range(len.split,len.orac),
     xlab="Number of path steps",ylab="Avg length",
     main=paste0("Split conformal + stepwise:\nAverage length"))
points(1:length(len.split),len.split,pch=20)
abline(h=len.orac,lty=2,col=2)
legend("topright",col=1:2,lty=1:2,
       legend=c("Split conformal","Oracle"))

# Plot test error
plot(1:length(err.split),err.split,type="l",
     ylim=range(err.split,err.orac),
     xlab="Number of path steps",ylab="Test error",
     main=paste0("Split conformal + stepwise:\nTest error"))
points(1:length(err.split),err.split,pch=20)
abline(h=err.orac,lty=2,col=2)
legend("topright",col=1:2,lty=1:2,
       legend=c("Split conformal","Oracle"))

####################

cat("Type return to continue ...\n")
tmp = readLines(n=1)

## Forward stepwise: cross-validate over 30 path steps

# Stepwise training and prediction functions
funs = lars.funs(type="stepwise",max.steps=30,cv=TRUE,cv.rule="1se")

# Split conformal inference 
out.split = conformal.pred.split(x, y, x0, alpha=0.1,
  train.fun=funs$train, predict.fun=funs$predict)

cov.split = mean(out.split$lo <= y0 & y0 <= out.split$up)
len.split = mean(out.split$up - out.split$lo)
err.split = mean((y0 - out.split$pred)^2)

# Compare to parametric intervals from oracle linear regression
lm.orac = lm(y~x[,1:s])
out.orac = predict(lm.orac, list(x=x0[,1:s]),
  interval="predict", level=0.9)

cov.orac = mean(out.orac[,"lwr"] <= y0 & y0 <= out.orac[,"upr"])
len.orac = mean(out.orac[,"upr"] - out.orac[,"lwr"])
err.orac = mean((y0 - out.orac[,"fit"])^2)

tab = matrix(c(cov.split,len.split,err.split,
  cov.orac,len.orac,err.orac),ncol=2)
colnames(tab) = c("Split conformal","Oracle")
rownames(tab) = c("Avg coverage","Avg length","Test error")
tab


