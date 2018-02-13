rm(list=ls())
gc()
source("functions/plots.R")

file.hig.hig <- 'output/compa-optim_p=100,s=10,n=200,eps=0.01.RData'
file.hig.med <- 'output/compa-optim_p=100,s=10,n=100,eps=0.01.RData'
file.hig.low <- 'output/compa-optim_p=100,s=10,n=50,eps=0.01.RData'
    
file.med.hig <- 'output/compa-optim_p=100,s=30,n=200,eps=0.01.RData'
file.med.med <- 'output/compa-optim_p=100,s=30,n=100,eps=0.01.RData'
file.med.low <- 'output/compa-optim_p=100,s=30,n=50,eps=0.01.RData'

pdf(file="output/timing_all.pdf", width=15, height=10)
plot.timer(file.med.low, bw=0.25, title=FALSE)
dev.off()

## Possible files...
file.other.low <- "output/quadra_vs_others_p=200,n=400.RData"
file.other.hig <- "output/quadra_vs_others_p=400,n=100.RData"

pdf(file="timing_others_high_dim.pdf")
plot.others(file.other.hig, title=FALSE)
dev.off()

pdf(file="timing_others_low_dim.pdf")
plot.others(file.other.low, title=FALSE)
dev.off()

## Possible files...
file.mine.hig <- "output/quadra_all_p=400,n=100.RData"

pdf(file="timing_quadra_all_high_dim.pdf")
plot.mines(file.mine.hig, title=FALSE)
dev.off()

## final plot for packages comparison

file.other.low <- "output/timing_others_p=100,n=40_nsim=30.RData"
file.other.med <- "output/timing_others_p=1000,n=200_nsim=30.RData"
file.other.hig <- "output/timing_others_p=10000,n=400_nsim=10.RData"

pdf(file="timing_others_low.pdf")
par(mar=c(2,2,0.1,0.1))
plot.others(file.other.low, title=FALSE, legend=FALSE, cex=1.5)
dev.off()

pdf(file="timing_others_med.pdf")
par(mar=c(2,2,0.1,0.1))
plot.others(file.other.med, title=FALSE, legend=FALSE, cex=1.5)
dev.off()

pdf(file="timing_others_hig.pdf")
par(mar=c(2,2,0.1,0.1))
plot.others(file.other.hig, title=FALSE, legend=FALSE, cex=1.5)
dev.off()

pdf(file="timing_others_legend.pdf")
par(mar=c(2,2,0.1,0.1))
plot(NA,xlim=c(-1,1),ylim=-c(-1,1), axes=FALSE, xlab="", ylab="")
legend("topleft", c("low correlation (0.1)","med correlation (0.4)","high correlation (0.8)"), seg.len=3,
       lty=c(1,1,1), col=c("blue","green","red"), lwd=2, bty="n", cex=2)
legend("bottomleft", c("glmnet (CD, active set)", "SPAMs (FISTA, no active set)", "SPAMs (homotopy/LARS)", "quadrupen (this paper)", "lars (homotopy/LARS)"),
       lty=c(2,3,-1,-1,-1), seg.len=3, lwd=2, pch = c(16, 17, 17, 15, 16), col="black", bty="n", cex=2)
dev.off()
