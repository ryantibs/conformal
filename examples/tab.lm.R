# Print linear model simulations, low-dimensional

dim = "lo"
method.nums = 1:4
metric.nums = c(1,2,5)
method.names = c("Conformal","Jackknife","Split","Parametric")
metric.names = c("Coverage", "Length","Time")
lets = c("A","B","C")
tun = rep(1,4) # Which tuning parameter value for each method?

for (j in 1:length(lets)) {
  name = paste0("lm.",dim,".sim",lets[j])
  x = readRDS(paste0("rds/",name,".rds"))

  tab = tab.se = matrix(0,5,length(method.nums)) 
  colnames(tab) = method.names
  for (i in 1:length(method.nums)) {
    num = method.nums[i]
    tab[,i] = c(x$ave.cov[[num]][tun[i]],x$ave.len[[num]][tun[i]],
         x$ave.err[[num]][tun[i]],x$ave.opt[[num]][tun[i]],
         x$ave.tim[[num]])
    tab.se[,i] = c(x$sd.cov[[num]][tun[i]],x$sd.len[[num]][tun[i]],
            x$sd.err[[num]][tun[i]],x$sd.opt[[num]][tun[i]],
            x$sd.tim[[num]])/sqrt(nrow(x$cov[[num]]))
  }
  tab = tab[metric.nums,]; tab.se = tab.se[metric.nums,]
  rownames(tab) = metric.names
  print.tex(tab, tab.se, file=paste0("tab/",name,".tex"))
  cat(paste0("tab/",name,".tex\n"))
}

#####
# Print linear model simulations, high-dimensional, no ridge

dim = "hi"
method.nums = c(1,2,4)
metric.nums = c(1,2,5)
method.names = c("Conformal","Jackknife","Parametric")
metric.names = c("Coverage", "Length","Time")
lets = c("A","B","C")
tun = rep(1,3) # Which tuning parameter value for each method?

for (j in 1:length(lets)) {
  name = paste0("lm.",dim,".sim",lets[j])
  x = readRDS(paste0("rds/",name,".rds"))

  tab = tab.se = matrix(0,5,length(method.nums)) 
  colnames(tab) = method.names
  for (i in 1:length(method.nums)) {
    num = method.nums[i]
    tab[,i] = c(x$ave.cov[[num]][tun[i]],x$ave.len[[num]][tun[i]],
         x$ave.err[[num]][tun[i]],x$ave.opt[[num]][tun[i]],
         x$ave.tim[[num]])
    tab.se[,i] = c(x$sd.cov[[num]][tun[i]],x$sd.len[[num]][tun[i]],
            x$sd.err[[num]][tun[i]],x$sd.opt[[num]][tun[i]],
            x$sd.tim[[num]])/sqrt(nrow(x$cov[[num]]))
  }
  tab = tab[metric.nums,]; tab.se = tab.se[metric.nums,]
  rownames(tab) = metric.names
  print.tex(tab, tab.se, file=paste0("tab/",name,".tex"))
  cat(paste0("tab/",name,".tex\n"))
}

#####
# Print linear model simulations, high-dimensional, ridge

dim = "hi"
method.nums = c(1,2,5,4)
metric.nums = c(1,2,3,5)
method.names = c("Conformal","Jackknife","Split+CV","Parametric")
metric.names = c("Coverage","Length","Test error","Time")
lets = c("A","B","C")
tun = c(2,2,1,2) # Which tuning parameter value for each method?

for (j in 1:length(lets)) {
  name = paste0("lm.",dim,".sim",lets[j])
  x = readRDS(paste0("rds/",name,".rds"))

  tab = tab.se = matrix(0,5,length(method.nums)) 
  colnames(tab) = method.names
  for (i in 1:length(method.nums)) {
    num = method.nums[i]
    tab[,i] = c(x$ave.cov[[num]][tun[i]],x$ave.len[[num]][tun[i]],
         x$ave.err[[num]][tun[i]],x$ave.opt[[num]][tun[i]],
         x$ave.tim[[num]])
    tab.se[,i] = c(x$sd.cov[[num]][tun[i]],x$sd.len[[num]][tun[i]],
            x$sd.err[[num]][tun[i]],x$sd.opt[[num]][tun[i]],
            x$sd.tim[[num]])/sqrt(nrow(x$cov[[num]]))
  }
  tab = tab[metric.nums,]; tab.se = tab.se[metric.nums,]
  rownames(tab) = metric.names
  print.tex(tab, tab.se, file=paste0("tab/ridge.",name,".tex"))
  cat(paste0("tab/",name,".tex\n"))
}
