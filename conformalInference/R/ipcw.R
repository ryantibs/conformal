#' Inverse Probability Censoring Weights
#' 
#' Computation of Inverse Probability Censoring Weights. The censoring survival 
#' function can either be estimated with Kaplan-Meier, Cox or Random Survival 
#' Forests.
#'
#' @param x Matrix of variables, of dimension (say) n x p.
#' @param t Vector of responses (observed times), of length (say) n.
#' @param d Vector of responses (censoring indicator), of length (say) n.
#' @param tau Horizon of time.
#' @param cens.model Model used to estimate the censoring survival function.
#'   This should be equal to "km", "rsf" or "cox" (respectively Kaplan-Meier, 
#'   Random Survival Forests or Cox). Default is "km".
#'
#' @return A vector of length n containing the censoring weights.
#'
#' @export ipcw

ipcw = function(t,d,x,tau,cens.model = "km"){
  # if(!is.null(x)){  x = as.matrix(x)  }
  if(!cens.model %in% c("km","cox","rsf")) stop("cens.model must be either km, 
                                                cox or rsf")
  
  if(cens.model == "km"){
    return(ipcw.km(t,d,x,tau))
  }else if(cens.model == "cox"){
    return(ipcw.cox(t,d,x,tau))
  }else if(cens.model == "rsf"){
    return(ipcw.rfsrc(t,d,x,tau))
  }
}

ipcw.km = function(t,d,x,tau){
  n = length(t)
  df = data.frame(t=t,d=1-d)
  fit = survfit(Surv(t, d) ~ 1, data=df)
  censfct = stepfun(fit$time,c(1,fit$surv))
  G.hat.tau = censfct(tau)
  w = as.double(lapply(1:n, function(s) {
    if (t[s] <= tau) {
      if (d[s] && censfct(t[s]) == 0) {
        warning("Error in weight calculation: Division by zero")
        return(NA)
      }
      return(d[s] / ifelse(d[s],censfct(t[s]),1))
    } else {
      if (G.hat.tau == 0) {
        warning("Error in weight calculation: Division by zero")
        return(NA)
      }
      return(1 / G.hat.tau)
    }
  }))
  return(w)
}

# ipcw.cox = function(t,d,x,tau){
#   n = length(t)
#   df = data.frame(t=t,d=1-d,x)
#   fit = survfit(coxph(Surv(t, d) ~ ., data=df), newdata=data.frame(x))
#   censfct = sapply(1:n,function(j) stepfun(fit$time,c(1,fit$surv[,j])))
#   w = as.double(lapply(1:n, function(s) {
#     if (t[s] <= tau) {
#       if (d[s] && censfct[[s]](t[s]) == 0) {
#         warning("Error in weight calculation: Division by zero")
#         return(NA)
#       }
#       return(d[s] / ifelse(d[s],censfct[[s]](t[s]),1))
#     } else {
#       if (censfct[[s]](tau) == 0) {
#         warning("Error in weight calculation: Division by zero")
#         return(NA)
#       }
#       return(1 / censfct[[s]](tau))
#     }
#   }))
#   return(w)
# }

ipcw.cox = function(t,d,x,tau){
  n = length(t)
  df = data.frame(t=t,d=1-d,x)
  coxfitC = coxph(Surv(t, d) ~ ., data=df)
  cumlambdaC = basehaz(coxfitC,centered=FALSE) #cumulative baseline hazard
  hazC=diff(c(0,cumlambdaC$hazard)) #baseline hazard
  timeC=cumlambdaC$time[hazC!=0] #times where the hazard is not null
  hazC=hazC[hazC!=0] #hazard for the times where the hazard is not null
  baseC_fun = stepfun(timeC,c(0,cumsum(hazC)))
  G.hat.t = exp(-baseC_fun(t)*exp(as.matrix(x)%*%coxfitC$coefficients))
  G.hat.tau = exp(-baseC_fun(tau)*exp(as.matrix(x)%*%coxfitC$coefficients))
  return(ifelse(t<=tau,d/G.hat.t,1/G.hat.tau))
}

ipcw.rfsrc = function(t,d,x,tau){
  n = length(t)
  df = data.frame(t=t,d=1-d,x)
  fit = rfsrc(Surv(t, d) ~ ., data=df)
  censfct = sapply(1:n,function(j) stepfun(fit$time.interest,c(1,fit$survival[j,])))
  w = as.double(lapply(1:n, function(s) {
    if (t[s] <= tau) {
      if (d[s] && censfct[[s]](t[s]) == 0) {
        warning("Error in weight calculation: Division by zero")
        return(NA)
      }
      return(d[s] / ifelse(d[s],censfct[[s]](t[s]),1))
    } else {
      if (censfct[[s]](tau) == 0) {
        warning("Error in weight calculation: Division by zero")
        return(NA)
      }
      return(1 / censfct[[s]](tau))
    }
  }))
  return(w)
}
