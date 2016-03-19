# Simulation setting B: now, nonlinear mean, and t2 errors
cat("========================================\n")
cat("                SETTING B               \n")
cat("========================================\n")
simB = sim.master(n, p, conformal.pred.funs, n0=n0, nrep=nrep, verb=TRUE,
  file=paste0(path,"simB.rds"), x.dist="normal", cor="none",
  mean.fun="additive", m=4, s=s, error.dist="t", df=2, sigma=sigma, snr=snr,
  sigma.type="const")
