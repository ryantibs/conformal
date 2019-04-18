## Statistics
load("rda/airfoil_sim.rda")
apply(cov.mat, 2, mean, na.rm=TRUE)
apply(len.mat, 2, median, na.rm=TRUE)
colMeans(!is.finite(len.mat))
apply(beta.mat1, 2, mean, na.rm=TRUE)
apply(beta.mat2, 2, mean, na.rm=TRUE)
mean(ssize)

## Histograms, coverage
xlim = range(cov.mat[,1:6], na.rm=TRUE)
bks = seq(xlim[1],xlim[2],length=60)
ylim = c()
for (i in 1:6) {
  ylim = range(c(ylim, hist(cov.mat[,i],breaks=bks,plot=FALSE)$count))
}

hh = 8; ww = 6
mar = c(4.25,4.25,1,1)
xlab = "Coverage"; ylab = "Frequency"

pdf(file="fig/airfoil_cov.pdf", height=hh, width=ww)
par(mar=mar, mfrow=c(3,1))
hist(cov.mat[,1], breaks=bks, col="red", xlim=xlim, ylim=ylim,
     main="", xlab=xlab, ylab=ylab)
hist(cov.mat[,2], breaks=bks, col=rgb(0,0,1,0.5), add=TRUE)
rug(mean(cov.mat[,1]), col="red", lwd=4)
rug(mean(cov.mat[,2]), col=rgb(0,0,1,0.5), lwd=4)
legend("topleft", col=c("red",rgb(0,0,1,0.5)), lwd=4, lty=1,
       legend=c("No covariate shift", "Covariate shift"))

hist(cov.mat[,3], breaks=bks, col="orange", xlim=xlim, ylim=ylim, 
     main="", xlab=xlab, ylab=ylab)
hist(cov.mat[,4], breaks=bks, col=rgb(0.5,0,0.5,0.5), add=TRUE)
rug(mean(cov.mat[,3]), col="orange", lwd=4)
rug(mean(cov.mat[,4]), col=rgb(0.5,0,0.5,0.5), lwd=4)
legend("topleft", col=c("orange",rgb(0.5,0,0.5,0.5)), lwd=4, lty=1,
       legend=c("Oracle weights", "No shift, fewer samples"))

hist(cov.mat[,5], breaks=bks, col="gray65", xlim=xlim, ylim=ylim, 
     main="", xlab=xlab, ylab=ylab)
hist(cov.mat[,6], breaks=bks, col=rgb(0,1,0,0.5), add=TRUE)
rug(mean(cov.mat[,5]), col="gray65", lwd=4)
rug(mean(cov.mat[,6]), col=rgb(0,1,0,0.5), lwd=4)
legend("topleft", col=c("gray65",rgb(0,1,0,0.5)), lwd=4, lty=1,
       legend=c("Logistic regression weights", "Random forest weights"))
graphics.off()

## Histograms, length
xlim = range(len.mat[,1:6], na.rm=TRUE, finite=TRUE)
bks = seq(xlim[1],xlim[2],length=75)
ylim = c()
for (i in 1:6) {
  ylim = range(c(ylim, hist(len.mat[,i],breaks=bks,plot=FALSE)$count))
}

hh = 8; ww = 6
mar = c(4.25,4.25,1,1)
main = ""
xlab = "Length"; ylab = "Frequency"
xlim = c(xlim[1],40) # Manually trim this, due to RF long right tail

pdf(file="fig/airfoil_len.pdf", height=hh, width=ww)
par(mar=mar, mfrow=c(3,1))
hist(len.mat[,1], breaks=bks, col="red", xlim=xlim, ylim=ylim,
     main="", xlab=xlab, ylab=ylab)
hist(len.mat[,2], breaks=bks, col=rgb(0,0,1,0.5), add=TRUE)
rug(median(len.mat[,1]), col="red", lwd=4)
rug(median(len.mat[,2]), col=rgb(0,0,1,0.5), lwd=4)
legend("topright", col=c("red",rgb(0,0,1,0.5)), lwd=4, lty=1,
       legend=c("No covariate shift", "Covariate shift"))

hist(len.mat[,3], breaks=bks, col="orange", xlim=xlim, ylim=ylim, 
     main="", xlab=xlab, ylab=ylab)
hist(len.mat[,4], breaks=bks, col=rgb(0.5,0,0.5,0.5), add=TRUE)
rug(median(len.mat[,3]), col="orange", lwd=4)
rug(median(len.mat[,4]), col=rgb(0.5,0,0.5,0.5), lwd=4)
legend("topright", col=c("orange",rgb(0.5,0,0.5,0.5)), lwd=4, lty=1,
       legend=c("Oracle weights", "No shift, fewer samples"))

hist(len.mat[,5], breaks=bks, col="gray65", xlim=xlim, ylim=ylim, 
     main="", xlab=xlab, ylab=ylab)
hist(len.mat[,6], breaks=bks, col=rgb(0,1,0,0.5), add=TRUE)
rug(median(len.mat[,5]), col="gray65", lwd=4)
rug(median(len.mat[,6]), col=rgb(0,1,0,0.5), lwd=4)
legend("topright", col=c("gray65",rgb(0,1,0,0.5)), lwd=4, lty=1,
       legend=c("Logistic regression weights", "Random forest weights"))
graphics.off()

## Histograms, adaptive?
xlim1 = range(cov.mat[,c(1,7:8)], na.rm=TRUE, finite=TRUE)
bks1 = seq(xlim1[1],xlim1[2],length=40)
ylim1 = c()
for (i in c(1,7:8)) {
  ylim1 = range(c(ylim1, hist(cov.mat[,i],breaks=bks1,plot=FALSE)$count))
}
xlab1 = "Coverage"; ylab1 = "Frequency"; cex = 0.9

xlim2 = range(len.mat[,c(1,7:8)], na.rm=TRUE, finite=TRUE)
bks2 = seq(xlim2[1],xlim2[2],length=70)
ylim2 = c()
for (i in c(1,7:8)) {
  ylim2 = range(c(ylim2, hist(len.mat[,i],breaks=bks2,plot=FALSE)$count))
}
xlab2 = "Length"; ylab2 = "Frequency"
xlim2 = c(xlim2[1],25) # Manually trim this, due to RF long right tail

hh = 6.5; ww = 6.5
mar = c(4.25,4.25,1,1)

pdf(file="fig/airfoil_ada.pdf", height=hh, width=ww)
par(mar=mar, mfrow=c(2,1))
hist(cov.mat[,1], breaks=bks1, col="red", xlim=xlim1, ylim=ylim1,
     main="", xlab=xlab1, ylab=ylab1, cex.lab=cex, cex.axis=cex)
hist(cov.mat[,7], breaks=bks1, col=rgb(0.65,0.65,0.65,0.5), add=TRUE)
hist(cov.mat[,8], breaks=bks1, col=rgb(0,1,0,0.5), add=TRUE)
rug(mean(cov.mat[,1]), col="red", lwd=4)
rug(mean(cov.mat[,7]), col="gray65", lwd=4)
rug(mean(cov.mat[,8]), col=rgb(0,1,0,0.5), lwd=4)
legend("topleft", col=c("red", rgb(0.65,0.65,0.65,0.5),rgb(0,1,0,0.5)),
       lwd=4, lty=1, legend=c("Unity weights", "Logistic regression weights",
                              "Random forest weights"), cex=cex)

hist(len.mat[,1], breaks=bks2, col="red", xlim=xlim2, ylim=ylim2,
     main="", xlab=xlab2, ylab=ylab2, cex.lab=cex, cex.axis=cex)
hist(len.mat[,7], breaks=bks2, col=rgb(0.65,0.65,0.65,0.5), add=TRUE)
hist(len.mat[,8], breaks=bks2, col=rgb(0,1,0,0.5), add=TRUE)
rug(median(len.mat[,1]), col="gray65", lwd=4)
rug(median(len.mat[,7]), col="gray65", lwd=4)
rug(median(len.mat[,8]), col=rgb(0,1,0,0.5), lwd=4)
graphics.off()

## Tilted densities
set.seed(0)
i = sample(N,n)
x = dat.x[i,]; y = dat.y[i]
x0 = dat.x[-i,]; y0 = dat.y[-i]
i0 = wsample(w(x0))
n00 = length(i0); x00 = x0[i0,]; y00 = y0[i0]

cex = 1.3
hh = 8; ww = 8.5
mar = c(4.25,4.25,1,1)

pdf(file="fig/airfoil_tilt.pdf", height=hh, width=ww)
par(mar=mar, mfcol=c(2,2))
plot(x0[,1], y0, main="", xlab="x1", ylab="y",
     cex.lab=cex, cex.axis=cex)

d1 = density(x0[,1]); d2 = density(x00[,1])
plot(d1, type="l", ylim=c(min(d1$y,d2$y),max(d1$y,d2$y)),
     main="", xlab="x1", ylab="Frequency",
     cex.lab=cex, cex.axis=cex)
lines(d2, col="blue")

plot(x0[,5], y0, main="", xlab="x5", ylab="y",
     cex.lab=cex, cex.axis=cex)

d1 = density(x0[,5]); d2 = density(x00[,5])
plot(d1, type="l", ylim=c(min(d1$y,d2$y),max(d1$y,d2$y)),
     main="", xlab="x5", ylab="Frequency",
     cex.lab=cex, cex.axis=cex)
lines(d2, col="blue")
graphics.off()
