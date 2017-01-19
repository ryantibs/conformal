# Simulation setting A: classical setup --- linear (fully supported) mean
# function; normal, homoskedastic errors; and normal, uncorrelated features
cat("========================================\n")
cat("                SETTING A               \n")
cat("========================================\n")
simA = sim.master(n, p, conformal.pred.funs, n0=n0, nrep=nrep, verb=TRUE,
  file=paste0(path,"simA.rds"), x.dist="normal", cor="none",
  mean.fun="linear", s=s, error.dist="normal", sigma=sigma, bval=bval,
  sigma.type="const")
