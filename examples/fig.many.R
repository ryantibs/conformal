# Plot all comparisons of conformal methods, across settings A through E,
# low- and high-dimensional

method.nums = c(1,2,5,3,4)
method.names = c("Lasso","Elastic net","Stepwise","Spam","Random forest")
cols = c(1,2,"purple",3,4)
lets = c("A","B","C","D","E")
w = 5
h = 5.5

for (dim in c("lo","hi")) {
for (j in 1:length(lets)) {
  name = paste0("many.",dim,".sim",lets[j])
  x = readRDS(paste0("rds/",name,".rds"))
  main = paste0(c("Coverage","Length","Test Error"),", Setting ",lets[j])
    
  plot(x, se=T, log="x", method.nums=method.nums, method.names=method.names,
       main=main, cols=cols, make.pdf=T, fig.dir="fig/", file.prefix=name,
       cex.main=1.5, w=w, h=h)

  cat(paste0("fig/",name,".pdf\n"))
}}

