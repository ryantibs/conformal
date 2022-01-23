## Lasso: use a fixed sequence of 100 lambda values

# Generate some example training data
set.seed(33)
n = 200; p = 100; s = 10
x = matrix(rnorm(n*p),n,p)
beta = c(rnorm(s),rep(0,p-s)) 
y = x %*% beta + rnorm(n)

# Generate some example test data
n0 = 100
x0 = matrix(rnorm(n0*p),n0,p)
y0 = x0 %*% beta + rnorm(n0)

# Grab a fixed lambda sequence from one call to glmnet, then
# define lasso training and prediction functions
out.gnet = glmnet::glmnet(x,y,nlambda=100,lambda.min.ratio=1e-4)
lambda = out.gnet$lambda
funs = lasso.funs(lambda=lambda)


# Split conformal inference 
out.multisplit = conformal.pred.msplit(x, y, x0, alpha=0.1,
                                       train.fun=funs$train, predict.fun=funs$predict, B = 40)

y0.mat = matrix(rep(y0,ncol(out.multisplit$lo)),nrow=n0)
cov.multisplit = colMeans(out.multisplit$lo <= y0.mat & y0.mat <= out.multisplit$up)
len.multisplit = colMeans(out.multisplit$up - out.multisplit$lo)

#### Compare to parametric intervals from oracle linear regression
lm.orac = lm(y~x[,1:s])
out.orac = predict(lm.orac,newdata=list(x=x0[,1:s]),
                   interval="predict", level=0.9)

cov.orac = mean(out.orac[,"lwr"] <= y0 & y0 <= out.orac[,"upr"])
len.orac = mean(out.orac[,"upr"] - out.orac[,"lwr"])

# Plot average coverage 
plot(log(lambda),cov.multisplit,type="o",pch=20,ylim=c(0,1),
     xlab="log(lambda)",ylab="Avg coverage",
     main=paste0("Multi Split + lasso (fixed lambda sequence):",
                 "\nAverage coverage"))
abline(h=cov.orac,lty=2,col=2)
legend("bottomleft",col=1:2,lty=1:2,
       legend=c("Multi Split","Oracle"))

# Plot average length
plot(log(lambda),len.multisplit,type="o",pch=20,
     ylim=range(len.multisplit,len.orac),
     xlab="log(lambda)",ylab="Avg length",
     main=paste0("Multi Split + lasso (fixed lambda sequence):",
                 "\nAverage length"))
abline(h=len.orac,lty=2,col=2)
legend("topleft",col=1:2,lty=1:2,
       legend=c("Multi Split","Oracle"))

#### Compare to split conformal intervals


out.split = conformal.pred.split(x, y, x0, alpha=0.1,
                                 train.fun=funs$train, predict.fun=funs$predict)

y0.mat = matrix(rep(y0,ncol(out.split$lo)),nrow=n0)
cov.split = colMeans(out.split$lo <= y0.mat & y0.mat <= out.split$up)
len.split = colMeans(out.split$up - out.split$lo)


# Plot average coverage 
plot(log(lambda),cov.multisplit,type="o",pch=20,ylim=c(0,1),
     xlab="log(lambda)",ylab="Avg coverage", col = "blue",
     main=paste0("Multi Split + lasso (fixed lambda sequence):",
                 "\nAverage coverage"))
points(log(lambda),cov.split, col="red", lty=2, pch=18)
legend("bottomleft",col=c("blue","red"),
       legend=c("Multi Split","Split Conformal"), lty = 1:2)

# Plot average length
plot(log(lambda),len.multisplit,type="o",pch=20,
     ylim=range(len.multisplit,len.split),
     xlab="log(lambda)",ylab="Avg length",
     main=paste0("Multi Split + lasso (fixed lambda sequence):",
                 "\nAverage length"))
points(log(lambda),len.split, col="red", lty=2, pch=18)
legend("bottomleft",col=c("blue","red"),
       legend=c("Multi Split","Split Conformal"), lty = 1:2)




####################

cat("Type return to continue ...\n")
tmp = readLines(n=1)

## Lasso: alternate way, use a dynamic sequence of 100 lambda values

# Lasso training and prediction functions
funs = lasso.funs(nlambda=100)

# multisplit conformal inference 
out.multisplit = conformal.pred.multisplit(x, y, x0, alpha=0.1,
                                       train.fun=funs$train, predict.fun=funs$predict,
                                       B=40)

y0.mat = matrix(rep(y0,ncol(out.multisplit$lo)),nrow=n0)
cov.multisplit = colMeans(out.multisplit$lo <= y0.mat & y0.mat <= out.multisplit$up)
len.multisplit = colMeans(out.multisplit$up - out.multisplit$lo)

# Compare to parametric intervals from oracle linear regression
lm.orac = lm(y~x[,1:s])
out.orac = predict(lm.orac,newdata=list(x=x0[,1:s]),
                   interval="predict", level=0.9)

cov.orac = mean(out.orac[,"lwr"] <= y0 & y0 <= out.orac[,"upr"])
len.orac = mean(out.orac[,"upr"] - out.orac[,"lwr"])

# Plot average coverage 
plot(log(length(cov.multisplit):1),cov.multisplit,type="o",pch=20,ylim=c(0,1),
     xlab="log(lambda rank) (i.e., log(1)=0 is smallest)",
     ylab="Avg coverage",
     main=paste0("Multi Split + lasso (dynamic lambda sequence):",
                 "\nAverage coverage"))
abline(h=cov.orac,lty=2,col=2)
legend("bottomleft",col=1:2,lty=1:2,
       legend=c("Multi Split","Oracle"))

# Plot average length
plot(log(length(len.multisplit):1),len.multisplit,type="o",pch=20,
     ylim=range(len.multisplit,len.orac),
     xlab="log(lambda rank) (i.e., log(1)=0 is smallest)",
     ylab="Avg length",
     main=paste0("Multi Split + lasso (dynamic lambda sequence):",
                 "\nAverage length"))
abline(h=len.orac,lty=2,col=2)
legend("topleft",col=1:2,lty=1:2,
       legend=c("Multi Split","Oracle"))



#### Compare to conformal split intervals
out.split = conformal.pred.split(x, y, x0, alpha=0.1,
                                 train.fun=funs$train, predict.fun=funs$predict)

y0.mat = matrix(rep(y0,ncol(out.split$lo)),nrow=n0)
cov.split = colMeans(out.split$lo <= y0.mat & y0.mat <= out.split$up)
len.split = colMeans(out.split$up - out.split$lo)
err.split = colMeans((y0.mat - out.split$pred)^2)

# Plot average coverage 
plot(log(length(cov.multisplit):1),cov.multisplit,type="o",pch=20,ylim=c(0,1),
     xlab="log(lambda rank) (i.e., log(1)=0 is smallest)",
     ylab="Avg coverage",
     main=paste0("Multi Split + lasso (dynamic lambda sequence):",
                 "\nAverage coverage"), col ="blue")
points(log(length(cov.split):1),cov.split, col="red", lty=2, pch=18)
legend("bottomleft",col=c("blue","red"),
       legend=c("Multi Split","Split Conformal"), lty = 1:2)

# Plot average length
plot(log(length(len.multisplit):1),len.multisplit,type="o",pch=20,
     ylim=range(len.multisplit,len.orac),
     xlab="log(lambda rank) (i.e., log(1)=0 is smallest)",
     ylab="Avg length",
     main=paste0("Multi Split + lasso (dynamic lambda sequence):",
                 "\nAverage length"), col="blue")
points(log(length(len.split):1),len.split, col="red", lty=2, pch=18)
legend("bottomleft",col=c("blue","red"),
       legend=c("Multi Split","Split Conformal"), lty = 1:2)

###########################################################################
# What if I varied the value of tau ?

funs =lasso.funs(cv=TRUE)

B=40

out.multisplit = conformal.pred.multisplit(x, y, x0, alpha=0.1,
                                           train.fun=funs$train, predict.fun=funs$predict,
                                           B=B, tau = 0)


out.multisplit2 = conformal.pred.multisplit(x, y, x0, alpha=0.1,
                                           train.fun=funs$train, predict.fun=funs$predict,
                                           B=B, tau = 1-(B+1)/(6*B))

avg_len = mean(out.multisplit$up - out.multisplit$lo)
avg_len

avg_len2 = mean(out.multisplit2$up - out.multisplit2$lo)
avg_len2

#As tau increases, the average interval length decreases