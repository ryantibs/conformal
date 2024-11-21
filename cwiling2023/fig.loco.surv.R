#############
# LIBRARIES #
#############


libraries_sources = function(){
  library(survival)
  library(FastPseudo)
  library(pseudo)
  library(randomForest)
  library(randomForestSRC)
  library(tidyverse)
  library(gridExtra)
  library(haven)
  library(devtools)
  library(knitr)
  library(Hmisc)
  library(parallel)
  library(conformalInference)
}
libraries_sources()

########################
# SIMULATION FUNCTIONS #
########################

simulate.A = function(n,seed=NULL,p=2,a=3,beta_tilde0 = c(5.5,2.5,2.5), 
                       beta0 = c(5.5,2.097,2.097,3.16), lambda=0.07,
                       C.indep.Z=T, k=12,nu=6,b=c(2,1)){
  if(!is.null(seed)){set.seed(seed)}
  x = matrix(rbinom(n*p,1,0.5), nrow=n)
  t = cbind(rep(1,n), x) %*% beta_tilde0 + runif(n,-a,a)
  if(C.indep.Z){c = rexp(n,rate=lambda)}
  else{c = k*(-log(1-runif(n))*exp(-as.double(x%*%b)))**(1/nu)}
  tobs = pmin(t,c)
  delta = (t <= c)*1
  mu_tau = cbind(rep(1,n),x[,1]*(1-x[,2]),x[,2]*(1-x[,1]),x[,1]*x[,2])%*% beta0
  return(list(t=t,c=c,tobs=tobs,delta=delta,x=x,mu_tau=mu_tau))
}


simulate.B = function(n,seed=NULL,a=5,k=2,nu=6,p=3,b=c(2,1,0),lambda=0.3){
  if(!is.null(seed)){set.seed(seed)}
  x = matrix(runif(n*p,-a,a),n,p)
  u = runif(n)
  t = k*(-log(1-u)*exp(-as.double(x%*%b)))**(1/nu)
  c = rexp(n,lambda)
  tobs = pmin(t,c)
  delta = (t <= c)*1
  return(list(t=t,c=c,tobs=tobs,delta=delta,x=x))
}


simulate.C = function(n,seed=NULL,a=5,k=2,nu=6,p.bin=6,p.cont=9,lambda=0.3){
  if(!is.null(seed)){set.seed(seed)}
  xbin = matrix(rbinom(n*p.bin,1,prob=0.4),n,p.bin)
  xcont = matrix(runif(n*p.cont,0,1),n,p.cont)
  x = cbind(xcont[,1],xbin[,1],xcont[,2],xbin[,2],xcont[,3],xbin[,3],xcont[,c(4,5)],xbin[,4],xcont[,6],xbin[,-seq(1,4)],xcont[,-seq(1,6)])
  y = x[,3] - 3*x[,5] + 2*x[,1]*x[,10] + 4*x[,2]*x[,7] + 3*x[,4]*x[,5] - 5*x[,6]*x[,10] + 3*x[,8]*x[,9] + x[,1]*x[,4] - 2*x[,6]*x[,9] - 4*x[,3]*x[,4] - x[,7]*x[,8]
  u = runif(n)
  t = k*(-log(1-u)*exp(-as.double(y)))**(1/nu)
  c = rexp(n,lambda)
  tobs = pmin(t,c)
  delta = (t <= c)*1
  return(list(t=t,c=c,tobs=tobs,delta=delta,x=x))
}


###############################
# TRAIN AND PREDICT FUNCTIONS #
###############################

integral = function(f,min=0,max){
  return(integrate(f,min,max,subdivisions=3000,rel.tol = 1e-5)$value)
}

train.fun = function(x,t,d,tau,pred.models=c("all"),tstar=NULL) {
  bincols = colMeans((x==1|x==0), na.rm=T) == 1
  for(i in 1:length(bincols)){
    # Turn the binary columns into factors
    if(bincols[i]){x[[i]] = as.factor(x[[i]])}
  }
  datafr = data.frame(x)
  
  pred.list = list()
  
  # Kaplan-Meier
  if("km" %in% pred.models || "all" %in% pred.models){
    pred.list$out.km = survfit(Surv(time=t,event=d)~1, data=datafr)
  }
  # Cox
  if("cox" %in% pred.models || "all" %in% pred.models){
    pred.list$out.cox = coxph(Surv(time=t,event=d)~., data = datafr)
  }
  # Random Survival Forests
  if("rsf" %in% pred.models || "all" %in% pred.models){
    pred.list$out.rsf = rfsrc(Surv(time=t,event=d)~.,data=data.frame(t=t,d=d,x),
                              importance=T, splitrule="bs.gradient", mtry=ncol(datafr))
  }
  # Pseudo-observations + Linear Model
  if("lm" %in% pred.models || "all" %in% pred.models){
    pobs = pseudoKM(t,as.logical(d),tau)$pseudoval
    pred.list$out.lm = lm(pobs~.,data=datafr)
  }
  # True event times + Linear Model with interactions (for scheme A only)
  if(("lm.oracle" %in% pred.models) && !is.null(tstar)){
    tstar.min.tau = pmin(tstar,tau)
    dat = data.frame(X01=(x[,1]==0 & x[,2]==1), 
                     X10=(x[,1]==1 & x[,2]==0),
                     X11=(x[,1]==1 & x[,2]==1))
    pred.list$out.lm.oracle = lm(tstar.min.tau~X01+X10+X11,data=dat)
  }
  # Cox with interactions (adapted to scheme C)
  if("cox.oracle" %in% pred.models){
    fml = paste("Surv(time=t,event=d)~(",paste0("X",1:ncol(datafr),collapse = "+"),")^2")
    pred.list$out.cox.oracle = coxph(as.formula(fml), data = datafr)
  }
  
  return(pred.list)
}

predict.fun = function(out,x0,tau,pred.models=c("all")) {
  predictions = numeric(0)
  
  # Kaplan-Meier
  if("km" %in% pred.models || "all" %in% pred.models){
    km.curve = stepfun(out$out.km$time,c(1,out$out.km$surv))
    predictions = c(predictions,rep(integral(km.curve,0,tau),nrow(x0)))
  }
  # Cox
  if("cox" %in% pred.models || "all" %in% pred.models){
    surv.cox = survfit(out$out.cox,newdata=data.frame(x0))
    pred.cox = numeric(nrow(x0))
    for(i in 1:nrow(x0)){
      cox.curve = stepfun(surv.cox$time,c(1,surv.cox$surv[,i]))
      pred.cox[i] = integral(cox.curve,0,tau)
    }
    predictions = c(predictions,pred.cox)
  }
  # Random Survival Forests
  if("rsf" %in% pred.models || "all" %in% pred.models){
    obj.pred.rsf = predict(out$out.rsf,data.frame(x0))
    vals_list = split(obj.pred.rsf$survival, row(obj.pred.rsf$survival))
    surv.rsf = purrr::map(vals_list, ~ stepfun(obj.pred.rsf$time.interest,c(1,.),f=1))
    pred.rsf = purrr::map_dbl(surv.rsf, ~integral(.,0,tau))
    predictions = c(predictions,as.double(pred.rsf))
  }
  # Pseudo-observations + Linear Model
  if("lm" %in% pred.models || "all" %in% pred.models){
    predictions = c(predictions,as.double(predict(out$out.lm,data.frame(x0))))
  }
  # True event times + Linear Model with interactions (for scheme A only)
  if("lm.oracle" %in% pred.models){
    dat = data.frame(X01=(x0[,1]==0 & x0[,2]==1), 
                     X10=(x0[,1]==1 & x0[,2]==0),
                     X11=(x0[,1]==1 & x0[,2]==1))
    predictions = c(predictions,as.double(predict(out$out.lm.oracle,dat)))
  }
  # Cox with interactions (adapted to scheme C)
  if("cox.oracle" %in% pred.models){
    surv.cox = survfit(out$out.cox.oracle,newdata=data.frame(x0))
    pred.cox = numeric(nrow(x0))
    for(i in 1:nrow(x0)){
      cox.curve = stepfun(surv.cox$time,c(1,surv.cox$surv[,i]))
      pred.cox[i] = integral(cox.curve,0,tau)
    }
    predictions = c(predictions,as.double(pred.cox))
  }
  
  return(predictions)
}

###################
# PARALLELISATION #
###################

parallelisation_init = function(){
  cl <<- makeCluster(getOption("cl.cores", detectCores()))
  invisible(clusterEvalQ(cl, {
    library(survival)
    library(FastPseudo)
    library(randomForest)
    library(randomForestSRC)
    library(tidyverse)
    library(gridExtra)
    library(haven)
    library(devtools)
    library(knitr)
    library(Hmisc)
    library(conformalInference)
  }))
  clusterExport(cl=cl,c("simulate.A","simulate.B","simulate.C",
                        "integral","train.fun","predict.fun"))
}


###########################################
# 1/WRSS: CONVERGENCE TO MSE(\Tilde{\mu}) #
###########################################

mse = function(n,n0,models,cens.model="km"){
  # Data simulation
  df = simulate(n+n0)
  df$w = ipcw(df$tobs,df$delta,df$x,tau,cens.model)
  #Training set
  i1 = sample(1:(n+n0),n,replace=F)
  df1 = list(t=df$t[i1],c=df$c[i1],tobs=df$tobs[i1],delta=df$delta[i1],
             x=df$x[i1,],w=df$w[i1])
  #Test set
  i0 = (1:(n+n0))[-i1]
  df0 = list(t=df$t[i0],c=df$c[i0],tobs=df$tobs[i0],delta=df$delta[i0],
             x=df$x[i0,],w=df$w[i0])
  # Train on full data set and fit on new data set
  out = train.fun(df1$x,df1$tobs,df1$delta,tau,models,df1$t) 
  fit = matrix(predict.fun(out,df0$x,tau,models),nrow=n0)
  WRSS = apply(df0$w*(pmin(df0$tobs,tau) - fit)**2,2,sum)/n0
  return(WRSS)
}

mse.pll = function(x,cens.model="km"){
  set.seed(x)
  return(purrr::cross_df(list(n=N)) %>%
           mutate(WRSS = map(n, ~ mse(.,.,models,cens.model))) )
}


##############################################
# 1.1/ CENSORING INDEPENDENT FROM COVARIATES #
##############################################

#Schema A1 (Linear model with two covariates; accessible RMST; C, Z independent)
simulate = simulate.A; tau = 8.8
models = c("lm","lm.oracle")
models.names = c("P.obs.+LM","Oracle")
m = length(models.names)
######

## Monte Carlo estimation of coefficients beta0
# df1 = simulate(5e7,seed=0)
# b00 = mean(pmin(df1$t[df1$x[,1]==0 & df1$x[,2]==0],tau))
# b10 = mean(pmin(df1$t[df1$x[,1]==1 & df1$x[,2]==0],tau)) - b00
# b01 = mean(pmin(df1$t[df1$x[,1]==0 & df1$x[,2]==1],tau)) - b00
# b11 = mean(pmin(df1$t[df1$x[,1]==1 & df1$x[,2]==1],tau)) - b00
# c(b00,b01,b10,b11)
##

## Monte Carlo estimation of the inseparability term
# inseparability.term = mean((pmin(df1$t,tau) - df1$mu_tau)**2)
inseparability.term = 1.5759
##

## Monte Carlo estimation of the imprecision term for the model P.obs.+LM
# train.fun.lm = function(x,t,d,tau) {
#   pobs= pseudoKM(t,as.logical(d),tau)$pseudoval
#   out.lm = lm(pobs~X1+X2,data=data.frame(x))
#   return(out.lm)
# }
# coeffs_A1 = sapply(1:20, function(i){
#   df = simulate(1e4,seed=i)
#   return(train.fun.lm(df$x,df$tobs,df$delta,tau)$coefficients)
# })
# beta_tilde_A1 = as.double(apply(coeffs_A1,1,median))
# imprecision.term.A1 = mean((df1$mu_tau - cbind(rep(1,1e7), df1$x) %*% beta_tilde_A1)**2)
imprecision.term.A1 = 0.067
##

## Illustration of Theorem 1
N = c(100,500,1000)
nb.repeat = 1000

parallelisation_init()
clusterExport(cl,c("simulate","tau","models","N","mse"))
mse.estimations.1 = parLapply(cl, X = 1:nb.repeat, fun = mse.pll)

mse.estimations.1 = purrr::cross_df(list(model = models.names, n=N,
                                         rep = 1:nb.repeat)) %>%
  mutate(wrss = unlist(data.table::rbindlist(mse.estimations.1)$WRSS))

mse.estimations.1$model = factor(mse.estimations.1$model, levels=models.names)
stopCluster(cl)

pdf(file="fig/WRSS.indep.pdf",w=4,h=2.5)
mse.estimations.1 %>% 
  ggplot(aes(x = as.factor(n), y = wrss)) +
  geom_boxplot(lwd=0.2,outlier.size = 0.1) +
  geom_hline(yintercept = inseparability.term, linetype = "dashed", 
             color = "red", lwd = 0.3)+
  geom_hline(data = mse.estimations.1 %>% filter(model == "P.obs.+LM"),
             aes(yintercept = inseparability.term + imprecision.term.A1),
             linetype="dashed", color = "blue", lwd = 0.3) +
  labs(x = 'Train and test size', y = 'WRSS') +
  facet_grid( .  ~ fct_rev(model)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
graphics.off()


############################################
# 1.2/ CENSORING DEPENDENT FROM COVARIATES #
############################################

#Schema A2 (Linear model with two covariates;accessible RMST;C dependent from Z)
simulate = function(n,seed=NULL){simulate.A(n,seed=seed,C.indep.Z = F)}; tau = 8.8
models = c("lm","lm.oracle")
models.names = c("P.obs.+LM","Oracle")
m = length(models.names)
######

## Monte Carlo estimation of the imprecision term for the model P.obs.+LM
# coeffs_A2 = sapply(1:20, function(i){
#   df = simulate(1e4,seed=i)
#   return(train.fun.lm(df$x,df$tobs,df$delta,tau)$coefficients)
# })
# beta_tilde_A2 = as.double(apply(coeffs_A2,1,median))
# imprecision.term.A2 = mean((df1$mu_tau - cbind(rep(1,1e7), df1$x) %*% beta_tilde_A2)**2)
imprecision.term.A2 = 0.0835
##

## Illustration of Theorem 1
N = c(100,500,1000)
nb.repeat = 1000

parallelisation_init()
clusterExport(cl,c("simulate","tau","models","N","mse","mse.pll"))
mse.estimations.km = parLapply(cl, X = 1:nb.repeat, 
                                fun = function(x){mse.pll(x,cens.model="km")})

mse.estimations.cox = parLapply(cl, X = 1:nb.repeat, 
                                 fun = function(x){mse.pll(x,cens.model="cox")})

mse.estimations.rsf = parLapply(cl, X = 1:nb.repeat, 
                                 fun = function(x){mse.pll(x,cens.model="rsf")})
stopCluster(cl)

mse.estimations.2 = purrr::cross_df(list(model = models.names, n=N,
                                         rep = 1:nb.repeat, 
                                         cens.model = c("KM","Cox","RSF"))) %>%
  mutate(wrss = c(unlist(data.table::rbindlist(mse.estimations.km)$WRSS),
                  unlist(data.table::rbindlist(mse.estimations.cox)$WRSS),
                  unlist(data.table::rbindlist(mse.estimations.rsf)$WRSS)))

mse.estimations.2$model = factor(mse.estimations.2$model, levels=models.names)
mse.estimations.2$cens.model = factor(mse.estimations.2$cens.model, levels=c("KM","Cox","RSF"))



pdf(file="fig/WRSS.dep.pdf",w=5,h=4)
mse.estimations.2 %>% 
  ggplot(aes(x = as.factor(n), y = wrss)) +
  scale_y_continuous(limits = c(1,3.5)) +
  geom_boxplot(lwd=0.2,outlier.size = 0.1) +
  geom_hline(yintercept = inseparability.term, linetype = "dashed", 
             color = "red", lwd = 0.3)+
  geom_hline(data = mse.estimations.2 %>% filter(model == "P.obs.+LM"),
             aes(yintercept = inseparability.term + imprecision.term.A2),
             linetype="dashed", color = "blue", lwd = 0.3) +
  labs(x = 'Train and test size', y = 'WRSS') +
  facet_grid( fct_rev(model)  ~ cens.model) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
graphics.off()


#####################
# 2/SPLIT ALGORITHM #
#####################

#Scheme B (Cox model, 3 covariates)
simulate = simulate.B; tau = 3.6;
models = c("km","cox","rsf","lm")
models.names = c("KM","Cox","RSF","P.obs.+LM")
m = length(models.names)
train.fun.2 = function(x,t,d,tau){train.fun(x,t,d,tau,pred.models=models)}
predict.fun.2 = function(out,x0,tau){predict.fun(out,x0,tau,pred.models=models)}
######

######################################
# 2.1/ ILLUSTRATION OF THE INTERVALS #
######################################

# Validation sample
n0 = 10
df0 = simulate(n0,seed=25)
df0$t.min.tau = pmin(df0$t,tau)

# Training sample
n = 4000
df = simulate(n,seed=40)

# Split algorithm
set.seed(0)
out.split = conformal.pred.split.surv(df$x, df$tobs, df$delta, tau, df0$x, 
                                      train.fun.2, predict.fun.2, cens.model = "km", rho=0.5, verbose=T)

res.split = data.frame(id = rep(LETTERS[1:n0],m), nb = rep(1:n0,m), 
                       t = rep(df0$t.min.tau,m), 
                       pred = as.double(out.split$pred),
                       lower = as.double(out.split$lo),
                       upper = as.double(out.split$up),
                       model = rep(models.names,each=n0))
res.split$model = factor(res.split$model, levels = models.names)


pdf(file="fig/split_all.pdf",w=6,h=4)
res.split %>%
  ggplot(aes(x=id,y=pred,shape=model,fill=model)) +  
  scale_shape_manual(values = c(23,22,21,24)) +
  geom_point(colour="black", fill="black", size=1.1,
             position=position_dodge(.5)) +
  geom_hline(yintercept=tau, linetype="dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper),position = position_dodge(.5),
                width=0.6,size=0.3) +
  geom_segment(aes(y = t, x = nb-0.35, yend = t, xend = nb+0.35),
               colour="red",size=0.45) +
  labs(x = "Validation set", 
       # title = "Predictions intervals",
       y = "Time to event", shape = "Regression model") +
  guides(shape = guide_legend(reverse = TRUE)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_flip()
graphics.off()




###########################
# 2.2/ EMPIRICAL COVERAGE #
###########################

cover = function(n,n0=N0){
  df0 = simulate(n0)
  df0$t.min.tau = pmin(df0$t,tau)
  df = simulate(n)
  out = conformal.pred.split.surv(df$x, df$tobs, df$delta, tau, df0$x, 
                                  train.fun.2, predict.fun.2, alpha = alpha, cens.model = "km", 
                                  rho=0.5, verbose=F)
  return(lapply(1:length(alpha),function(a)
    apply((df0$t.min.tau>=out$lo[,,a] & df0$t.min.tau<=out$up[,,a]),2,mean)))
}

cover.pll = function(x){
  set.seed(x)
  return(purrr::cross_df(list(n=N)) %>% mutate(coverage = map(n, ~ cover(.))) )
}

N = c(100,500,1000)
N0 = 500 
alpha = c(0.2,0.1,0.05)
nb.repeat = 200

parallelisation_init()
clusterExport(cl,c("simulate","tau","models","train.fun.2","predict.fun.2",
                   "N","N0","alpha","cover"))
coverage = parLapply(cl, X = 1:nb.repeat, fun = cover.pll)

coverage = purrr::cross_df(list(model = models.names, alpha=alpha, n=N, rep = 1:nb.repeat)) %>%
  mutate(coverage = unlist(data.table::rbindlist(coverage)$coverage))

coverage$model = factor(coverage$model, levels = models.names)

stopCluster(cl)


### Plot ###

#define custom color scale
myColors = c("#00364e","#006db5","#008bff")
names(myColors) = levels(coverage$n)
custom_colors = scale_colour_manual(values = myColors)

pdf(file="fig/coverage_fig.pdf",w=7,h=3)
coverage %>% 
  ggplot(aes(x = as.factor(1-alpha), y = coverage)) +
  coord_cartesian(ylim = c(0.585, 1)) +
  geom_boxplot(aes(color=as.factor(n)),lwd=0.3,outlier.size = 0.3) +
  geom_hline(yintercept=0.8, linetype="dashed", color = "grey",lwd=0.35) +
  geom_hline(yintercept=0.9, linetype="dashed", color = "grey",lwd=0.35) +
  geom_hline(yintercept=0.95, linetype="dashed", color = "grey",lwd=0.35) +
  labs(x = 'Confidence level', y = 'Empirical coverage', color='Training size') +
  custom_colors +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_grid( . ~ as.factor(model))
graphics.off()

### Tab ###

split.tab = coverage %>% group_by(alpha,n,model) %>% 
  dplyr::summarise(mean = round(mean(coverage),3),
                   median = round(median(coverage),3),
                   Q1 = round(quantile(coverage,0.25),3),
                   Q3 = round(quantile(coverage,0.75),3)) 

pdf("fig/coverage_tab.pdf",h=nrow(split.tab)/3.5,w=6)
grid.table(split.tab)
dev.off()


##########################
# 3/ VARIABLE IMPORTANCE #
##########################

#Scheme B (Cox model, 3 covariates)
simulate = simulate.B; tau = 3.6; p=3
models = c("km","cox","rsf","lm")
models.names = c("KM","Cox","RSF","P.obs.+LM")
m = length(models.names)
train.fun.3 = function(x,t,d,tau){train.fun(x,t,d,tau,pred.models=models)}
predict.fun.3 = function(out,x0,tau){predict.fun(out,x0,tau,pred.models=models)}
######

vimpG.pvalues = function(n){
  D2 = simulate(n)
  out.loco = loco.surv(rbind(D1$x,D2$x), c(D1$tobs,D2$tobs), c(D1$delta,D2$delta), 
                       tau, train.fun.3, predict.fun.3, split=1:n, vars=0,
                       rho=0.5, bonf.correct=FALSE, seed=NULL, verbose=FALSE)
  return(do.call(rbind,map(out.loco$inf.sign, ~ .[,"P-value"])))
}

vimpG.pvalues.pll = function(x){
  set.seed(x)
  return(vimpG.pvalues(n))
}

######

n=500
D1 = simulate(n,seed=0)


#########################
# 3.1/ ESTIMATION OF pk #
#########################

# parallelisation_init()
# set.seed(10)
# n2=1e5
# D2 = simulate(n2)
# mu = train.fun.3(D1$x,D1$tobs,D1$delta,tau)
# pred = matrix(predict.fun.3(mu,D2$x,tau),nrow=n2,ncol=m)
# res = abs(pmin(D2$t,tau)-pred)
# clusterExport(cl,c("tau","models","train.fun.3","predict.fun.3",
#                    "D2","D1","n2","m","res"))
# pk = parSapply(cl, X = 1:p, FUN = function(k){
#   set.seed(k)
#   mu.k = train.fun.3(D1$x[,-k],D1$tobs,D1$delta,tau)
#   pred.k = matrix(predict.fun.3(mu.k,D2$x[,-k],tau),nrow=n2,ncol=m)
#   delta.k = abs(pmin(D2$t,tau)-pred.k) - res
#   return(apply(delta.k,2,FUN=function(s){mean((s>=0)[D2$t<=tau])}))
# })
# colnames(pk) = paste0("Z",1:3)
# rownames(pk) = models.names
# stopCluster(cl)

pk = matrix(c(1.0000000, 1.0000000, 1.0000000, 
              0.8710390, 0.7945491, 0.4876769, 
              0.8154177, 0.7148376, 0.4360612, 
              0.8436275, 0.6984591, 0.4578992),4,3,byrow=T)
colnames(pk) = paste0("Z",1:3)
rownames(pk) = models.names

sink(file = "fig/pk.txt", type = "output")
print(round(pk,3))
sink()


##################################################
# 3.2/ VARIABLE IMPORTANCE MEASURES ILLUSTRATION #
##################################################

set.seed(3)
D2 = simulate(n)
out.D1 = rmst.pred(rbind(D1$x,D2$x), c(D1$tobs,D2$tobs), c(D1$delta,D2$delta), 
                   tau, train.fun.3, predict.fun.3, split=1:n, varsL=0, varsG=0,
                   rho=0.5, bonf.correct=FALSE, seed=NULL, verbose=TRUE,
                   error=F, roo=F,vimpL=T, vimpG=T)

### Local variable importance ###

pdf(file="fig/vimpLD1.pdf",w=7.5,h=6.5)
par(mfcol=c(p,m))
plot(out.D1, model.names = models.names, elements = c("vimpL"))
graphics.off()

### Global variable importance ###

pdf(file="fig/vimpGD1.pdf",w=7,h=2.25)
par(mfrow=c(1,m))
plot(out.D1, model.names = models.names, elements = c("vimpG"))
graphics.off()


###########################
# 3.3/ CALIBRATION SLOPES #
###########################
nb.repeat = 1000

parallelisation_init()
clusterExport(cl,c("simulate","tau","models","train.fun.3","predict.fun.3",
                   "n","D1","vimpG.pvalues"))
pvalues = parLapply(cl, X = 1:nb.repeat, fun = vimpG.pvalues.pll)

pvalues = array(unlist(pvalues), dim = c(m,p,nb.repeat))

stopCluster(cl)

pdf(file="fig/calibration.pdf",w=7.2,h=4.8)
par(mfcol=c(p,m))
for(i in 1:m){
  for(j in 1:p){
    hist(pvalues[i,j,],probability = T, main = paste0(models.names[i],", variable ",j),
         xlab = "P-value", ylab = "Density",xlim=c(0,1),breaks=seq(0,1,by=0.05))
  }
}
graphics.off()

# Rejection rate at level 5%
pvalues.reject = data.frame(t(apply(pvalues < 0.05,MARGIN = c(1,2),mean)))
colnames(pvalues.reject) = models.names

sink(file = "fig/reject_rate.txt", type = "output")
print(pvalues.reject)
sink()

##########################
# 4/ RMST.PRED FRAMEWORK #
##########################

integral = function(f,min=0,max){
  return(integrate(f,min,max,subdivisions=3000,rel.tol = 1e-4)$value)
}

#################
# 4.1/ SCHEME B #
#################

#Scheme B (Cox model, 3 covariates)
simulate = simulate.B; tau = 3.6; p=3
models = c("km","cox","rsf","lm")
models.names = c("KM","Cox","RSF","P.obs.+LM")
m = length(models.names)
train.fun.4 = function(x,t,d,tau){train.fun(x,t,d,tau,pred.models=models)}
predict.fun.4 = function(out,x0,tau){predict.fun(out,x0,tau,pred.models=models)}
######

n = 4000
df = simulate(n,seed=40)

### Analysis framework ###

V = 10
out.fmw.B = rmst.pred(df$x, df$tobs, df$delta, tau, train.fun.4, predict.fun.4, 
                      n.folds=V, cens.model="km", alpha=0.1, varsL=0, varsG=O, 
                      seed=60, verbose=TRUE, rho=0.5)

### MSE ###

pdf(file="fig/mseB.pdf",w=3.5,h=3.5)
par(mfrow=c(1,1))
plot(out.fmw.B, model.names = models.names, elements = c("mse"))
graphics.off()

### Local variable importance ###

pdf(file="fig/vimpLB.pdf",w=7.5,h=6.5)
par(mfcol=c(p,m))
plot(out.fmw.B, model.names = models.names, elements = c("vimpL"))
graphics.off()

### Global variable importance ###

## Test ##

sink(file = "fig/testB.txt", type = "output")
print(out.fmw.B, model.names = models.names)
sink()

## Display ##

pdf(file="fig/vimpGB.pdf",w=7,h=2.25)
par(mfrow=c(1,m))
plot(out.fmw.B, model.names = models.names, elements = c("vimpG"))
graphics.off()


#################
# 4.2/ SCHEME C #
#################

#Scheme C (Cox model, 15 covariates)
simulate = simulate.C; tau = 2.8; p = 15
models = c("km","cox","rsf","lm","cox.oracle")
models.names = c("KM","Cox","RSF","P.obs.+LM","Oracle")
m = length(models.names)
train.fun.4.2 = function(x,t,d,tau){train.fun(x,t,d,tau,pred.models=models)}
predict.fun.4.2 = function(out,x0,tau){predict.fun(out,x0,tau,pred.models=models)}
######

n = 1000
df = simulate(n,seed=40)

### Analysis framework ###

V = 10
out.fmw.C = rmst.pred(df$x, df$tobs, df$delta, tau, train.fun.4.2, predict.fun.4.2, 
                      n.folds=V, cens.model="km", alpha=0.1, varsG=0, 
                      seed=0, verbose=TRUE, rho=0.5, vimpL = F)

### MSE ###

pdf(file="fig/mseC.pdf",w=3.5,h=3.5)
par(mfrow=c(1,1))
plot(out.fmw.C, model.names = models.names, elements = c("mse"))
graphics.off()

### Global variable importance ###

## Test ##

sink(file = "fig/testC.txt", type = "output")
print(out.fmw.C, model.names = models.names)
sink()

## Display ##

pdf(file="fig/vimpGC.pdf",w=9,h=2.25)
par(mfrow=c(1,m))
plot(out.fmw.C, model.names = models.names, elements = c("vimpG"))
graphics.off()

## Re-computation of global vimp with another seed
out.fmw.C2 = rmst.pred(df$x, df$tobs, df$delta, tau, train.fun.4.2, predict.fun.4.2, 
                       n.folds=V, cens.model="km", alpha=0.1, varsG=0, 
                       seed=10, verbose=TRUE, rho=0.5, vimpL = F, error = F, roo=F)

pdf(file="fig/vimpGC2.pdf",w=9,h=2.25)
par(mfrow=c(1,m))
plot(out.fmw.C2, model.names = models.names, elements = c("vimpG"))
graphics.off()


######################
# 5/ MULTIPLE SPLITS #
######################

#Scheme C (Cox model, 15 covariates)
simulate = simulate.C; tau = 2.8; p = 15
models = c("km","cox","rsf","lm","cox.oracle")
models.names = c("KM","Cox","RSF","P.obs.+LM","Oracle")
m = length(models.names)
train.fun.5 = function(x,t,d,tau){train.fun(x,t,d,tau,pred.models=models)}
predict.fun.5 = function(out,x0,tau){predict.fun(out,x0,tau,pred.models=models)}
######

vars.active = 1:15

multiple.splits.pll = function(X){
  set.seed(X)
  out = loco.surv(df$x, df$tobs, df$delta, tau, train.fun.5, predict.fun.5, 
                  vars=vars.active,rho=0.5)
  return(sapply(1:length(out$inf.sign), 
                FUN = function(l){out$inf.sign[[l]][,1]}))
}

n = 1000
df = simulate(n,seed=40)
N = 50

parallelisation_init()
clusterExport(cl,c("simulate","models","tau","train.fun.5","predict.fun.5",
                   "df","vars.active"))
multiple.C = parLapply(cl, X = 1:N, fun = multiple.splits.pll)

multiple.C = array(unlist(multiple.C), dim = c(length(vars.active),m,N))

stopCluster(cl)

# Romano & DiCiccio 2019 procedure 2
# twice the median p-value or twice the average p-value serves as a valid 
# p-value for the overall test.

signif.code = function(v) {
  code = rep("",length(v))
  code[v < 0.1] = "."
  code[v < 0.05] = "*"
  code[v < 0.01] = "**"
  code[v < 0.001] = "***"
  return(code)
}

pval.twice.med = round(apply(multiple.C,MARGIN=2,
                             function(x){apply(x,1,function(y){pmin(2*median(y),1)})}),3)

tab.pval.twice.med = paste0("X",vars.active)
for(l in 1:m){
  tab = round(pval.twice.med[,l],digits=3)
  code = signif.code(pval.twice.med[,l])
  tab = cbind(tab,code)
  colnames(tab) = c("P-value",models.names[l])
  tab.pval.twice.med = cbind(tab.pval.twice.med,tab)
}
colnames(tab.pval.twice.med)[1] = "Var"

sink(file = "fig/multi_testC2.txt", type = "output")
print(tab.pval.twice.med,quote=F)
print("Significance codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1",quote=F)
sink()


###########################
# 6/ REAL DATA : BRCANCER #
###########################

# The outcome of interest is the time to relapse, not the time to death
brcancer0 = read_dta(url("https://www.stata-press.com/data/r16/brcancer.dta"))
brcancer = as.data.frame(as.matrix(brcancer0))
brcancer = brcancer %>% select(-id, -x4, -x4b, -x5e)
brcancer$x2 = brcancer$x2 - 1

# Events
tobs.brcancer = brcancer$rectime
delta.brcancer = brcancer$censrec
tau.brcancer = quantile(tobs.brcancer,0.9,names=F)
x.brcancer = brcancer %>% select(-rectime,-censrec)
x.brcancer$hormon = factor(x.brcancer$hormon)
x.brcancer$x2 = factor(x.brcancer$x2)
x.brcancer$x4a = factor(x.brcancer$x4a)
p = dim(x.brcancer)[2]

## Impact of covariates on censoring
cens_cox = coxph(Surv(tobs.brcancer,1-delta.brcancer)~.,data=x.brcancer)
cens_cox
cens_rsf = rfsrc(Surv(tobs.brcancer,cens.brcancer)~.,data=data.frame(tobs.brcancer,cens.brcancer=1-delta.brcancer,x.brcancer),
                 importance=T,seed=-1)
cens_rsf$importance
##

#####
models = c("km","cox","rsf","lm")
models.names = c("KM","Cox","RSF","P.obs.+LM")
m = length(models.names)
train.fun.6 = function(x,t,d,tau){train.fun(x,t,d,tau,pred.models=models)}
predict.fun.6 = function(out,x0,tau){predict.fun(out,x0,tau,pred.models=models)}
######

out.brcancer = 
  rmst.pred(x.brcancer, tobs.brcancer, delta.brcancer, tau.brcancer, 
            train.fun.6, predict.fun.6, n.folds=20, cens.model="km", 
            alpha=0.1, varsL=c(2,4,5,6,7), varsG=0, verbose=TRUE, seed=20,
            error=T,roo=T,vimpL=T,vimpG=T,rho=0.5)

### MSE ###

pdf(file="fig/mse_brcancer.pdf",w=3.5,h=3.5)
par(mfrow=c(1,1))
plot(out.brcancer, model.names = models.names, elements = c("mse"))
graphics.off()


### Local variable importance ###

pdf(file="fig/vimpL_brcancer.pdf",w=9,h=6)
par(mfrow=c(m,5))
plot(out.brcancer, model.names = models.names, elements = c("vimpL"))
graphics.off()

### Global variable importance ###

## Test ##

sink(file = "fig/test_brcancer.txt", type = "output")
print(out.brcancer, model.names = models.names, var.names = colnames(x.brcancer))
sink()

## Display ##

pdf(file="fig/vimpG_brcancer.pdf",w=12,h=3)
par(mfrow=c(1,m))
plot(out.brcancer, model.names = models.names, elements = c("vimpG"))
graphics.off()


### Multiple splitting ###

multiple.splits.br.pll = function(X){
  set.seed(X)
  out = loco.surv(x.brcancer, tobs.brcancer, delta.brcancer, tau.brcancer, 
                  train.fun.6, predict.fun.6, vars=0,rho=0.5)
  return(sapply(1:length(out$inf.sign), 
                FUN = function(l){out$inf.sign[[l]][,1]}))
}

N = 40

parallelisation_init()
clusterExport(cl,c("simulate","models","train.fun.6","predict.fun.6",
                   "x.brcancer","tobs.brcancer","delta.brcancer","tau.brcancer"))

multiple.brcancer = parLapply(cl, X = 1:N, fun = multiple.splits.br.pll)

multiple.brcancer = array(unlist(multiple.brcancer), dim = c(p,m,N))

stopCluster(cl)



pval.twice.med.br = round(apply(multiple.brcancer,MARGIN=2,
                                function(x){apply(x,1,function(y){pmin(2*median(y),1)})}),3)

tab.pval.twice.med.br = colnames(x.brcancer)
for(l in 1:m){
  tab = round(pval.twice.med.br[,l],digits=3)
  code = signif.code(pval.twice.med.br[,l])
  tab = cbind(tab,code)
  colnames(tab) = c("P-value",models.names[l])
  tab.pval.twice.med.br = cbind(tab.pval.twice.med.br,tab)
}
colnames(tab.pval.twice.med.br)[1] = "Var"

sink(file = "fig/multi_test_brcancer.txt", type = "output")
print(tab.pval.twice.med.br,quote=F)
print("Significance codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1",quote=F)
sink()

