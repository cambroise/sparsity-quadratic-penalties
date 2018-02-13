require(ggplot2)
plot.timer <- function(files, main="Timing comparison", title=TRUE, bw=0.25) {

  for (file in files) {
    load(file)
    sparsity <- factor(switch(strsplit(file, ",")[[1]][2],
                              "s=10"  = 'high sparsity',
                              "s=30" = 'medium sparsity'))
    sample <- factor(switch(strsplit(file, ",")[[1]][3],
                            "n=200"  = 'p = n/2',
                            "n=100" = 'p=n',
                            "n=50" = 'p = 2*n'))

    xy <- as.matrix(expand.grid(lambdas, gammas))
    
    prox.all <- c()
    prox.all     <- data.frame(rbind(cbind(xy, log10(c(res.prox.small$timer/res.quad.small$timer))),
                                     cbind(xy, log10(c(res.prox.medium$timer/res.quad.medium$timer))),
                                     cbind(xy, log10(c(res.prox.large$timer/res.quad.large$timer)))))
    prox.all <- cbind(prox.all, factor(rep(c("small corr. (0.1)", "medium corr. (0.4)", "large corr. (0.8)"), each=nrow(xy))))
    names(prox.all) <- c("x", "y", "z", "cor")
    
    path.all <- c()
    path.all     <- data.frame(rbind(cbind(xy, log10(c(res.path.small$timer/res.quad.small$timer))),
                                     cbind(xy, log10(c(res.path.medium$timer/res.quad.medium$timer))),
                                     cbind(xy, log10(c(res.path.large$timer/res.quad.large$timer)))))
    path.all <- cbind(path.all, factor(rep(c("small corr. (0.1)", "medium corr. (0.4)", "large corr. (0.8)"), each=nrow(xy))))
    names(path.all) <- c("x", "y", "z", "cor")
    
    res.all <- rbind(prox.all, path.all)
    res.all <- cbind(res.all, factor(rep(c("proximal (fista)", "coordinate descent"), c(nrow(prox.all),nrow(path.all)))))
    names(res.all) <- c("x", "y", "z", "correlation", "algo")
    
    d <- ggplot(data=res.all, aes(x=log10(y), y=log10(x), z=z)) 
    d <- d + geom_tile(aes(fill=z)) + stat_contour(size=0.2, binwidth=bw)
    if (title){
      d <- d + opts(title=paste(main,"with",sparsity,"and",sample)) + labs(x="lambda (log-scale)",y="gamma (log-scale)")
    } else {
      d <- d + labs(x="",y="")
    }
    d <- d + facet_grid(algo ~ correlation)  + scale_fill_continuous(name="# times\n   faster", breaks= c(0,log10(3),log10(10),log10(30),log10(100),log10(300)), labels=c("1","3","10","30","100","300"))
    print(d)
  }
}

plot.others <- function(file, main="timing of quadratic solver vs other packages on a full solution path",
                        cex = 1, y.range=NA, legend=TRUE, title=TRUE) {

  load(file)

  time.glmn.low <- apply(time.glmn.low,1, median, na.rm=TRUE)
  prec.glmn.low <- apply(prec.glmn.low,1, median, na.rm=TRUE)
 
  time.glmn.med <- apply(time.glmn.med,1,median, na.rm=TRUE)
  prec.glmn.med <- apply(prec.glmn.med,1,median, na.rm=TRUE)

  time.glmn.hig <- apply(time.glmn.high,1,median, na.rm=TRUE)
  prec.glmn.hig <- apply(prec.glmn.high,1,median, na.rm=TRUE)
  
  time.craf <- median(c(time.craf.low, time.craf.med, time.craf.high), na.rm=TRUE)
  prec.craf <- median(c(prec.craf.low, prec.craf.med, prec.craf.high), na.rm=TRUE)

  time.spam <- median(c(time.spam.low, time.spam.med, time.spam.high), na.rm=TRUE)
  prec.spam <- median(c(prec.spam.low, prec.spam.med, prec.spam.high), na.rm=TRUE)

  time.spamf.low <- apply(time.spamf.low,1, median, na.rm=TRUE)
  prec.spamf.low <- apply(prec.spamf.low,1, median, na.rm=TRUE)
 
  time.spamf.med <- apply(time.spamf.med,1,median, na.rm=TRUE)
  prec.spamf.med <- apply(prec.spamf.med,1,median, na.rm=TRUE)

  time.spamf.hig <- apply(time.spamf.high,1,median, na.rm=TRUE)
  prec.spamf.hig <- apply(prec.spamf.high,1,median, na.rm=TRUE)
  
  time.lars <- median(c(time.lars.low, time.lars.med, time.lars.high), na.rm=TRUE)

  all.time <- c(time.glmn.low,time.glmn.med,time.glmn.hig, time.spam,
                time.spamf.low,time.spamf.med,time.spamf.hig, time.craf)
  all.prec <- c(prec.glmn.low,prec.glmn.med,prec.glmn.hig, prec.spam,
                prec.spamf.low,prec.spamf.med,prec.spamf.hig, prec.craf)
  
  rx <- log10(range(all.time,na.rm=TRUE))
  if (is.na(y.range)) {
    ry <- log10(range(all.prec,na.rm=TRUE))
  } else {
    ry <- y.range
  }
  
  plot(log10(time.glmn.low),log10(prec.glmn.low), col="blue",
       xlab="", ylab="", ylim=ry, xlim=rx, type="b", lty=2, pch=16, cex=cex)
  if (title) {
    title(main = main, xlab="CPU time in sec (log10)", ylab="distance to optimum (log10)")
  }
  lines(log10(time.glmn.med),log10(prec.glmn.med), col="green",type="b", pch=16, lty=2, cex=cex)
  lines(log10(time.glmn.hig),log10(prec.glmn.hig), col="red",type="b", pch=16, lty=2, cex=cex)

  lines(log10(time.spamf.low),log10(prec.spamf.low), col="blue",type="b", pch=17, lty=3, cex=cex)
  lines(log10(time.spamf.med),log10(prec.spamf.med), col="green",type="b", pch=17, lty=3, cex=cex)
  lines(log10(time.spamf.hig),log10(prec.spamf.hig), col="red",type="b", pch=17, lty=3, cex=cex)
  
  segments(min(rx)-1,log10(prec.craf),log10(time.craf),log10(prec.craf), lty=3, cex=cex)
  segments(log10(time.craf),min(ry)-1,log10(time.craf),log10(prec.craf), lty=3, cex=cex)
  points(log10(time.craf),log10(prec.craf), pch=15, cex=cex)
  segments(log10(time.lars),min(ry)-1,log10(time.lars),log10(prec.craf), lty=3, cex=cex)
  points(log10(time.lars),log10(prec.craf), pch=16, cex=cex)

  segments(min(rx)-1,log10(prec.spam),log10(time.spam),log10(prec.spam), lty=3, cex=cex)
  segments(log10(time.spam),min(ry)-1,log10(time.spam),log10(prec.spam), lty=3, cex=cex)
  points(log10(time.spam),log10(prec.spam), pch=17, cex=cex)
  
  if (legend) {
    legend("topright", c("low corr.","med corr.","high corr."), seg.len=3,
           lty=c(1,1,1), col=c("blue","green","red"), lwd=2, bty="n")

    legend("bottomright", c("glmnet (CD, active set)", "SPAMs (FISTA, no active set)", "SPAMs (homotopy/LARS)", "quadrupen (this paper)", "lars (homotopy/LARS)"),
           lty=c(2,3,-1,-1,-1), seg.len=3, lwd=2, pch = c(16, 17, 17, 15, 16), col="black", bty="n")
  }
}

plot.mines <- function(file, main="timing of quadratic solver vs other methods on a full solution path", legend=TRUE, title=TRUE) {

  load(file)

  time.path.low <- apply(time.path.low,1, median, na.rm=TRUE)
  prec.path.low <- apply(prec.path.low,1, median, na.rm=TRUE)
 
  time.path.med <- apply(time.path.med,1,median, na.rm=TRUE)
  prec.path.med <- apply(prec.path.med,1,median, na.rm=TRUE)

  time.path.hig <- apply(time.path.high,1,median, na.rm=TRUE)
  prec.path.hig <- apply(prec.path.high,1,median, na.rm=TRUE)

  time.prox.low <- apply(time.prox.low,1, median, na.rm=TRUE)
  prec.prox.low <- apply(prec.prox.low,1, median, na.rm=TRUE)
 
  time.prox.med <- apply(time.prox.med,1,median, na.rm=TRUE)
  prec.prox.med <- apply(prec.prox.med,1,median, na.rm=TRUE)

  time.prox.hig <- apply(time.prox.high,1,median, na.rm=TRUE)
  prec.prox.hig <- apply(prec.prox.high,1,median, na.rm=TRUE)
  
  time.craf <- median(c(time.craf.low, time.craf.med, time.craf.high), na.rm=TRUE)
  prec.craf <- median(c(prec.craf.low, prec.craf.med, prec.craf.high), na.rm=TRUE)
  
  time.lars <- median(c(time.lars.low, time.lars.med, time.lars.high), na.rm=TRUE)

  all.time <- c(time.path.low,time.path.med,time.path.hig,
                time.prox.low,time.prox.med,time.prox.hig, time.craf)
  all.prec <- c(prec.path.low,prec.path.med,prec.path.hig,
                prec.prox.low,prec.prox.med,prec.prox.hig, prec.craf)
  
  rx <- log10(range(all.time,na.rm=TRUE))
  ry <- log10(range(all.prec,na.rm=TRUE))

  plot(log10(time.path.low),log10(prec.path.low), col="blue",
       xlab="", ylab="", ylim=ry, xlim=rx, type="b", lty=2, pch=16)
  if (title) {
    title(main = main, xlab="CPU time in sec (log10)", ylab="distance to optimum (log10)")
  }
  lines(log10(time.path.med),log10(prec.path.med), col="green",type="b", pch=16, lty=2)
  lines(log10(time.path.hig),log10(prec.path.hig), col="red",type="b", pch=16, lty=2)

  lines(log10(time.prox.low),log10(prec.prox.low), col="blue",type="b", pch=17, lty=3)
  lines(log10(time.prox.med),log10(prec.prox.med), col="green",type="b", pch=17, lty=3)
  lines(log10(time.prox.hig),log10(prec.prox.hig), col="red",type="b", pch=17, lty=3)
  
  segments(min(rx)-1,log10(prec.craf),log10(time.craf),log10(prec.craf), lty=3)
  segments(log10(time.craf),min(ry)-1,log10(time.craf),log10(prec.craf), lty=3)
  points(log10(time.craf),log10(prec.craf), pch=15)
  segments(log10(time.lars),min(ry)-1,log10(time.lars),log10(prec.craf), lty=3)
  points(log10(time.lars),log10(prec.craf), pch=16)

  if (legend) {
    legend("topright", c("low corr.","med corr.","high corr."), seg.len=3,
           lty=c(1,1,1), col=c("blue","green","red"), lwd=2, bty="n")

    legend("bottomright", c("quadrupen  (CD, active set)", "quadrupen  (FISTA, active set)", "quadrupen (homotopy/LARS)", "lars (homotopy/LARS)"),
           lty=c(2,3,-1,-1), seg.len=3, lwd=2, pch = c(16, 17, 15, 16), col="black", bty="n")
  }
}

multiplot <- function(..., plotlist=NULL, cols) {
    require(grid)

    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)

    # Make the panel
    plotCols = cols                          # Number of columns of plots
    plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols

    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
    vplayout <- function(x, y)
        viewport(layout.pos.row = x, layout.pos.col = y)

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
        curRow = ceiling(i/plotCols)
        curCol = (i-1) %% plotCols + 1
        print(plots[[i]], vp = vplayout(curRow, curCol ))
    }

}
