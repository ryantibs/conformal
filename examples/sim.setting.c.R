# Simulation setting C: now, nonlinear mean, heteroskedastic t2 errors, and
# auto-correlated mix of x variables
cat("========================================\n")
cat("                SETTING C               \n")
cat("========================================\n")
simC = sim.master(n, p, conformal.pred.funs, n0=n0, nrep=nrep, verb=TRUE,
  file=paste0(path,"simC.rds"), x.dist="mix", cor="auto", k=5,
  mean.fun="linear", s=s, error.dist="t", df=2, sigma=sigma, snr=snr,
  sigma.type="var")
