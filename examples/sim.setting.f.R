# Simulation setting F: back to nonlinear mean, heteroskedastic t2 errors, now
# random-correlated mix of x variables, rotated mean, weak sparsity
cat("========================================\n")
cat("                SETTING F               \n")
cat("========================================\n")
simF = sim.master(n, p, conformal.pred.funs, n0=n0, nrep=nrep, verb=TRUE,
  file=paste0(path,"simF.rds"), x.dist="mix", cor="rand", k=5,
  mean.fun="rotated", m=4, s=s, sparsity="weak", error.dist="t", df=2,
  sigma=sigma, bval=bval, sigma.type="var")
