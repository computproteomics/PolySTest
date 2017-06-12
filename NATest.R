library(matrixStats)
library(qvalue)
MissValPDistr <- function(NumReps, PercNA) {
  p <- PercNA
  d <- NumReps
  D <- rep(0,d+1)
  # terms of binomial distribution
  binTerms <- NULL
  for (i in 0:d) {
    binTerms <- append(binTerms, choose(d,i)*p^i*(1-p)^(d-i))
  }
  for (i in 0:d) {
    for (j in 0:(d-i)) {
      D[i+1] = D[i+1] + binTerms[j+1]*binTerms[j+i+1]
    }
  }
  for (i in 1:d) {
    D[i+1] <- 2*D[i+1]
  }
  D
}

NumReps <- 10
A  <- matrix(rnorm(10000*2*NumReps),ncol=NumReps*2)
A[1:50,1:NumReps] <- A[1:50,1:NumReps]+1
A[51:100,1:NumReps] <- A[51:100,1:NumReps]-1
AA <- A

qs <- quantile(A,probs=seq(0,1,0.01))
pvals <- statis <-  matrix(NA,nrow(A),ncol=length(qs))
for (q in qs)
{
  A[A<q] <- NA
  NAPDistr <- MissValPDistr(NumReps,sum(is.na(A))/(nrow(A)*2*NumReps))
  statis[,which(q==qs)] <- (rowSums(!is.na(A[,1:NumReps])) - rowSums(!is.na(A[,(NumReps+1):(2*NumReps)])))
  pvals[,which(q==qs)] <- NAPDistr[abs(statis[,which(q==qs)])+1]
  # hist(pNAvalues,100)
  
}
par(mfrow=c(1,1))
hist(pvals,100)
## Implement something like mean over the lowest values
hist(rowMins(pvals),100)
hist(pvals[60,],100)

pvalues <- rowMins(pvals)*(NumReps+1)
plot(ecdf(pvalues),log="xy",xlim=c(1e-3,1),ylim=c(1e-3,1))
abline(0,1)

sum(pvalues>1)/length(pvalues)
pvalues[pvalues>1] <- 1

qvalues <- qvalue(pvalues)$qvalue
sumqs <- NULL
for (i in seq(0,1,0.001))
  sumqs <- append(sumqs,sum(qvalues<i)/length(qvalues))
plot(seq(0,1,0.001), sumqs,type="l",log="xy")
abline(0,1)
which(qvalues<0.05)


par(mfrow=c(2,4))
lowest <- which.min(rowMins(pvals))
highest <- which.max(rowMins(pvals))
plot(qs,pvals[lowest,],ylab="p-values",xlab="threshold")
plot(AA[lowest,],col=rep(1:2,each=NumReps),ylab="Intensities")
hist(pvals[lowest,],50, main="p-value distribution (quantiles)")
plot(qs,statis[lowest,])
plot(qs,pvals[highest,],ylab="p-values")
plot(AA[highest,],col=rep(1:2,each=NumReps),ylab="Intensities")
hist(pvals[highest,],50, main="p-value distribution (quantiles)")
plot(qs,statis[highest,])

## testing circlize plots
library(circlize)
nfeat <- 10
NumCond <- 5
data <- matrix(rnorm(200),ncol=5*(NumCond-1),nrow=nfeat)


cols <- rainbow(nfeat)

par(mar=rep(0,4))
circos.clear()
circos.par(cell.padding=c(0,0,0,0),canvas.xlim=c(-0.5,0.5),canvas.ylim=c(-1.5,1.5),
           track.margin=c(0,0.01),start.degree=90,gap.degree=4)
circos.initialize(1:(NumCond-1),xlim=c(1,nfeat))

for (t in 1:5) {
  tsign <- data[,seq(t,ncol(data),NumCond-1)]>0
  circos.trackPlotRegion(ylim=c(-3,2),track.height=1/7, bg.border="#777777", 
                         panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    #circos.text(x=mean(xlim), y=1.7,
    #            labels=name, facing = dd, cex=0.6,  adj = aa),
    xdiff <- (xlim[2]-xlim[1])/nfeat
    if(t == 1) {
    circos.text(mean(xlim), max(ylim)+5, name, facing = "inside", niceFacing = TRUE)
    circos.axis("top", labels = rownames(data),major.at=seq(1.5,nfeat-0.5,length=nfeat),minor.ticks=0,labels.cex = 0.4,labels.facing = "reverse.clockwise")
    }
    for (j in 1:nfeat) {
      if (tsign[j,i])
        circos.rect(xleft=xlim[1]+(j-1)*xdiff, ybottom=ylim[1],
                    xright=xlim[2]-(nfeat-j)*xdiff, ytop=ylim[2],
                    col = cols[j], border=cols[j]
        )
    }})
}
fccols <- redblue(1001)
circos.trackPlotRegion(ylim=c(-3,2),track.height=1/7, bg.border=NA, panel.fun = function(x, y) {
  name = get.cell.meta.data("sector.index")
  i = get.cell.meta.data("sector.numeric.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  #circos.text(x=mean(xlim), y=1.7,
  #            labels=name, facing = dd, cex=0.6,  adj = aa),
  xdiff <- (xlim[2]-xlim[1])/nfeat
  for (j in 1:nfeat) {
    circos.rect(xleft=xlim[1]+(j-1)*xdiff, ybottom=ylim[1],
                  xright=xlim[2]-(nfeat-j)*xdiff, ytop=ylim[2],
                  col = fccols[(data[j,i]-min(data))*500], border=0
      )
  }})
text(0,0,"Fold\nchange",cex=0.5)
