# Simulation setting D: back to classical setup, but without sparsity
cat("========================================\n")
cat("                SETTING D               \n")
cat("========================================\n")
simD = sim.master(n, p, conformal.pred.funs, n0=n0, nrep=nrep, verb=TRUE,
  file=paste0(path,"simD.rds"), x.dist="normal", cor="none",
  mean.fun="linear", s=s2, error.dist="normal", sigma=sigma, bval=bval,
  sigma.type="const")
