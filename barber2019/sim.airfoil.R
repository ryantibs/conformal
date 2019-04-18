## Load the airfoil data
dat = read.table("airfoil.txt")
dim(dat)
colnames(dat) = c("Frequency",
                  "Angle",
                  "Chord",
                  "Velocity",
                  "Suction",
                  "Sound")

dat.x = as.matrix(dat[,1:5])
dat.y = as.numeric(dat[,6])
dat.x[,1] = log(dat.x[,1]) # Log transform
dat.x[,5] = log(dat.x[,5]) # Log transform
N = nrow(dat.x); p = ncol(dat.x)

## Exponential tilting functions
w = function(x) {
  exp(x[,c(1,5)] %*% c(-1,1))
}

wsample = function(wts, frac=0.25) {
  n = length(wts)
  i = c()
  while(length(i) <= n*frac) {
    i = c(i, which(runif(n) <= wts/max(wts)))
  }
  return(i)
}

## Run the simulation
library(randomForest)
library("conformalInference")

set.seed(0)
R = 5000
cov.mat = len.mat = matrix(0,R,8)
beta.mat1 = beta.mat2 = matrix(0,R,p+1)
ssize = numeric(R)
n = round(N/2)
my.lm.funs = lm.funs()

for (r in 1:R) {
  if (r %% 20 == 0) cat(r,".. ")
  i = sample(N,n)
  x = dat.x[i,]; y = dat.y[i]
  x0 = dat.x[-i,]; y0 = dat.y[-i]

  # No tilting
  out1 = conformal.pred.split(x, y, x0, 
                              train.fun=my.lm.funs$train,
                              predict.fun=my.lm.funs$predict)

  cov.mat[r,1] = mean(out1$lo <= y0 & y0 <= out1$up)
  len.mat[r,1] = median(out1$up - out1$lo)
  
  # Tilting
  i0 = wsample(w(x0))
  x00 = x0[i0,]; y00 = y0[i0]
  i1 = out1$split; n1 = length(out1$split)
  ssize[r] = sum(w(x[-i1,]))^2/sum(w(x[-i1,])^2)
  
  # Unweighted
  out2 = conformal.pred.split(x, y, x00,
                              train.fun=my.lm.funs$train,
                              predict.fun=my.lm.funs$predict)

  cov.mat[r,2] = mean(out2$lo <= y00 & y00 <= out2$up)
  len.mat[r,2] = median(out2$up - out2$lo)
  
  # Weighted, oracle  
  out3 = conformal.pred.split(x, y, x00, w=w(rbind(x,x00)),
                              train.fun=my.lm.funs$train,
                              predict.fun=my.lm.funs$predict)

  cov.mat[r,3] = mean(out3$lo <= y00 & y00 <= out3$up)
  len.mat[r,3] = median(out3$up - out3$lo)

  # No tilting, smaller sample size
  i2 = (1:n)[-i1]
  i3 = c(i1, sample(i2, floor(ssize[r])))
  out4 = conformal.pred.split(x[i3,], y[i3], x0, split=1:n1,
                              train.fun=my.lm.funs$train,
                              predict.fun=my.lm.funs$predict)

  cov.mat[r,4] = mean(out4$lo <= y0 & y0 <= out4$up)
  len.mat[r,4] = median(out4$up - out4$lo)
  
  # Weighted, estimated with logistic regression
  zy = as.factor(c(rep(0,n),rep(1,nrow(x00))))
  zx = rbind(x,x00)
  obj.glm = glm(zy ~ zx, family="binomial")
  prob.glm = predict(obj.glm, type="response")
  wts.glm = prob.glm / (1-prob.glm)
  beta.mat1[r,] = coef(obj.glm)
  
  out5 = conformal.pred.split(x, y, x00, w=wts.glm,
                              train.fun=my.lm.funs$train,
                              predict.fun=my.lm.funs$predict)

  cov.mat[r,5] = mean(out5$lo <= y00 & y00 <= out5$up)
  len.mat[r,5] = median(out5$up - out5$lo)

  # Weighted, estimated with random forest
  obj.rf = randomForest(zx, zy)
  prob.rf = predict(obj.rf, type="prob")[,2]
  prob.rf = pmax(pmin(prob.rf, 0.99), 0.01)
  wts.rf = prob.rf / (1-prob.rf)
  
  out6 = conformal.pred.split(x, y, x00, w=wts.rf,
                              train.fun=my.lm.funs$train,
                              predict.fun=my.lm.funs$predict)

  cov.mat[r,6] = mean(out6$lo <= y00 & y00 <= out6$up)
  len.mat[r,6] = median(out6$up - out6$lo)

  # No tilting, weighted, estimated with logistic
  zy = as.factor(c(rep(0,n),rep(1,nrow(x0))))
  zx = rbind(x,x0)
  obj.glm = glm(zy ~ zx, family="binomial")
  prob.glm = predict(obj.glm, type="response")
  wts.glm = prob.glm / (1-prob.glm)
  beta.mat2[r,] = coef(obj.glm)
  
  out7 = conformal.pred.split(x, y, x0, w=wts.glm,
                              train.fun=my.lm.funs$train,
                              predict.fun=my.lm.funs$predict)

  cov.mat[r,7] = mean(out7$lo <= y0 & y0 <= out7$up)
  len.mat[r,7] = median(out7$up - out7$lo)

  # Weighted, estimated with random forest
  obj.rf = randomForest(zx, zy)
  prob.rf = predict(obj.rf, type="prob")[,2]
  prob.rf = pmax(pmin(prob.rf, 0.99), 0.01)
  wts.rf = prob.rf / (1-prob.rf)
  
  out8 = conformal.pred.split(x, y, x0, w=wts.rf,
                              train.fun=my.lm.funs$train,
                              predict.fun=my.lm.funs$predict)

  cov.mat[r,8] = mean(out8$lo <= y0 & y0 <= out8$up)
  len.mat[r,8] = median(out8$up - out8$lo)
}
cat("\n")

save(list=ls(), file="rda/airfoil_sim.rda")
