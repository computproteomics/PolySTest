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
