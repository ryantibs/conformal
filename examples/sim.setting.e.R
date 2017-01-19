# Simulation setting E: back to classical setup, but with correlated variables,
# and omitted variable bias
cat("========================================\n")
cat("                SETTING E               \n")
cat("========================================\n")
simE = sim.master(n, p, conformal.pred.funs, n0=n0, nrep=nrep, verb=TRUE,
  file=paste0(path,"simE.rds"), x.dist="normal", cor="pair", rho=0.8,
  mean.fun="linear", s=s, error.dist="normal", sigma=sigma, bval=bval,
  sigma.type="const", omitted.vars=TRUE)
