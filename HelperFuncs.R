# Permutations: min 10 replicates, if not available, then generate from entire data


library(parallel)
library(limma)
source("rankprodbounds.R")

NumThreads <- 4
NTests <- 100 # for permutation tests
NumPermCols <- 7 # minimum number of columns for permutation tests (to ensure sufficient combinations)

StatsForPermutTest <- function(Data, Paired) {
  if (Paired) {
    # NumReps <- ncol(Data)
    Stats <- apply(Data, 1, function(x) mean(x,na.rm=T)/sd(x,na.rm=T))
    Stats <- abs(Stats) *sqrt(rowSums(!is.na(Data)))
  } else {
    NumReps <- ncol(Data)*.5
    NumRealReps <- rowSums(!is.na(Data))*0.5
    # pooled standard deviation
    Stats <- apply(Data, 1, function(x) (mean(x[1:NumReps],na.rm=T)-mean(x[(NumReps+1):(NumReps*2)],na.rm=T))/
                     (sqrt((var(x[1:NumReps],na.rm=T)+var(x[(NumReps+1):(NumReps*2)],na.rm=T)))))
    Stats <- abs(Stats)/((NumRealReps-1)*(2*NumRealReps))*((2*NumRealReps-2)*(NumRealReps*NumRealReps)) 
  }
  return(Stats)
}

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

RPStats <- function(tRPMAData,NumReps) {
  RPMADown_pvalues <- NULL
  NumElements<-rowSums(!is.na(tRPMAData))
  Rank <- NULL
  RP.own <- 0
  
  RPMADown_pvalues <- parallel::mclapply(unique(NumElements),function (d) {
    tRPMADown_pvalues <- NULL
    RPMAData<-tRPMAData[NumElements==d,]
    if(d>0 && length(as.matrix(RPMAData))>ncol(tRPMAData)) {
      RP.own<-0
      Rank<-NULL
      
      RankNAs<-0
      for (r in 1:NumReps) {
        Rank[[r]]<-rank(RPMAData[,r],na.last="keep")#/(sum(!is.na(RPMAData[,r]))+1)
        names(Rank[[r]]) <-rownames(RPMAData)
        Rank[[r]][is.na(Rank[[r]])]<-1
        RP.own<-RP.own+log(Rank[[r]])
        RankNAs<-RankNAs+sum(Rank[[r]]>1)
      }
      RP.own<-exp(RP.own)
      tNumReps <- d
      tNumFeat <- length(RP.own)
      # print(range(RP.own))
      # print(tNumFeat)
      # print(RP.own)
      RP.own <- round(RP.own)
      tRPout <- rankprodbounds(RP.own, tNumFeat, tNumReps, "geometric")
      names(tRPout) <- names(RP.own)
      tRPMADown_pvalues <- c(tRPMADown_pvalues,tRPout)
      tRPMADown_pvalues
    }
  },mc.cores=NumThreads)
  # print(head(unlist(RPMADown_pvalues)))
  return(unlist(RPMADown_pvalues))
}


Paired <- function(MAData,NumCond,NumReps) {
  ##########################################################
  MAReps<-rep(1:NumCond,NumReps)
  
  Rank <- NULL
  
  ##limma with ratios
  design<-plvalues<-NULL
  for (c in (1:(NumCond))) {
    design<-cbind(design,as.numeric(MAReps==c))
  }
  lm.fittedMA <- lmFit(MAData,design)
  lm.bayesMA<-eBayes(lm.fittedMA)
  topTable(lm.bayesMA)
  plvalues <- lm.bayesMA$p.value
  qlvalues <- matrix(NA,nrow=nrow(plvalues),ncol=ncol(plvalues),dimnames=dimnames(plvalues))
  # qvalue correction
  for (i in 1:ncol(plvalues)) {
    tqs <- qvalue(na.omit(plvalues[,i]))$qvalues
    qlvalues[names(tqs),i] <- tqs
  }
  
  ## rank products + t-test + permutation test
  ptvalues<-NULL
  pRPvalues<-matrix(NA,ncol=NumCond,nrow=nrow(MAData),dimnames=list(rows = rownames(MAData), cols=paste("RP p-values",1:NumCond)))
  pPermutvalues<-matrix(NA,ncol=NumCond,nrow=nrow(MAData),dimnames=list(rows = rownames(MAData), cols=paste("Permutation p-values",1:NumCond)))
  for (vs in 1:NumCond) {
    if (!is.null(getDefaultReactiveDomain()))
      setProgress(0.1+0.3/NumCond*vs, detail = paste("tests for comparison",vs,"of",NumCond))
    
    tMAData<-MAData[,MAReps==vs]
    ptMAvalues<-NULL
    ## MA t-test_pvalues
    for (pep in 1:(dim(tMAData)[1])) {
      if(sum(!is.na(tMAData[pep,]))>1) {
        ptMAvalues<-append(ptMAvalues,t.test(as.vector(tMAData[pep,]))$p.value)
      } else {
        ptMAvalues <- append(ptMAvalues,NA)
      }
    }
    names(ptMAvalues)<-rownames(tMAData)
    ptvalues <- cbind(ptvalues,ptMAvalues)
    
    ## rank products
    tRPMAData<-MAData[,MAReps==vs]
    #Up
    RPMAUp_pvalues <- RPStats(tRPMAData,NumReps)
    #Down
    RPMADown_pvalues <- RPStats(-tRPMAData,NumReps)
    ttt <- rowMins(cbind(RPMAUp_pvalues,RPMADown_pvalues),na.rm=T)*2
    ttt[ttt>1] <- 1
    pRPvalues[names(RPMAUp_pvalues),vs] <- ttt
    
    ## Permutation tests: add columns from randomized full set to reach min. NumPermCols replicates
    # randomizing also sign to avoid tendencies to one or the other side
    if (ncol(tMAData)<NumPermCols) {
      AddDat <- matrix(sample(as.vector(tMAData),(NumPermCols-ncol(tMAData))*nrow(tMAData),replace=T),nrow=nrow(tMAData))
      PermMAData <- cbind(tMAData, AddDat)
    } else {
      PermMAData <- tMAData
    }
    RealStats <- StatsForPermutTest(tMAData,Paired=T)
    PermutOut <- parallel::mclapply(1:NTests,function (x) {indat <- apply(PermMAData,1,function(y) sample(y,NumReps)*sample(c(1,-1),NumReps,replace=T))
    StatsForPermutTest(t(indat),T)
    },mc.cores=NumThreads)
    PermutOut <- matrix(unlist(PermutOut),nrow=nrow(tMAData))
    # print(cbind(RealStats,PermutOut)[1,])
    pPermutvalues[,vs] <- apply(cbind(RealStats,PermutOut), 1 , function(x) ifelse(is.na(x[1]),NA,(1+sum(x[1] < x[-1],na.rm=T))/(sum(!is.na(x)))))
    # print(PermMAData[1,])
    # print(head(pPermutvalues[,vs]))
  }
  lratios <- NULL
  qRPvalues <- qtvalues <- qPermutvalues <- matrix(NA,nrow=nrow(MAData),ncol=NumCond,dimnames=list(rows=rownames(MAData), cols=1:NumCond))
  for (i in 1:NumCond) {
    tqs <- qvalue(na.omit(ptvalues[,i]))$qvalues
    qtvalues[names(tqs),i] <- tqs
    print(range(pPermutvalues[,i]))
    tqs <- qvalue(na.omit(pPermutvalues[,i]))$qvalues
    qPermutvalues[names(tqs),i] <- tqs
    # print(range(na.omit(pRPvalues[,i])))
    tqs <- qvalue(na.omit(pRPvalues[,i]),lambda=seq(0.05,max(na.omit(pRPvalues[,i]))-0.05,0.05))$qvalues
    qRPvalues[names(tqs),i] <- tqs
    lratios <- cbind(lratios, rowMeans(MAData[,MAReps==i],na.rm=T))
  }
  
  
  
  return(list(lratios=lratios,ptvalues=ptvalues, plvalues=plvalues, pRPvalues=pRPvalues,pPermutvalues=pPermutvalues,
              qtvalues=qtvalues, qlvalues=qlvalues, qRPvalues=qRPvalues, qPermutvalues=qPermutvalues,Sds=sqrt(lm.bayesMA$s2.post)))
}

Unpaired <- function(Data,NumCond,NumReps) {
  ##########################################################
  # significance analysis
  Reps <- rep(1:NumCond,NumReps)
  
  ## limma
  design <- model.matrix(~0+factor(Reps-1))
  colnames(design)<-paste("i",c(1:NumCond),sep="")
  contrasts<-NULL
  First <- 1
  for (i in (1:NumCond)[-First]) contrasts<-append(contrasts,paste(colnames(design)[i],"-",colnames(design)[First],sep=""))
  contrast.matrix<-makeContrasts(contrasts=contrasts,levels=design)
  print(dim(Data))
  lm.fitted <- lmFit(Data,design)
  
  lm.contr <- contrasts.fit(lm.fitted,contrast.matrix)
  lm.bayes<-eBayes(lm.contr)
  topTable(lm.bayes)
  plvalues <- lm.bayes$p.value
  qlvalues <- matrix(NA,nrow=nrow(plvalues),ncol=ncol(plvalues),dimnames=dimnames(plvalues))
  # qvalue correction
  for (i in 1:ncol(plvalues)) {
    tqs <- qvalue(na.omit(plvalues[,i]))$qvalues
    qlvalues[names(tqs),i] <- tqs
  }
  
  ## Permutation tests: try individual an full data set???
  
  
  ## rank products + t-test
  ptvalues<-NULL
  pRPvalues<-matrix(NA,ncol=NumCond-1,nrow=nrow(Data),dimnames=list(rows = rownames(Data), cols=paste("RP p-values",1:(NumCond-1))))
  pPermutvalues<-matrix(NA,ncol=NumCond-1,nrow=nrow(Data),dimnames=list(rows = rownames(Data), cols=paste("Permutation p-values",1:(NumCond-1))))
  for (vs in 2:NumCond) {
    if (!is.null(getDefaultReactiveDomain()))
      setProgress(0.1+0.3/(NumCond-1)*vs, detail = paste("tests for comparison",vs-1,"of",NumCond-1))
    tData<-Data[,Reps==vs]
    trefData <- Data[,Reps==1]
    tptvalues<-NULL
    ## MA t-test_pvalues
    for (pep in 1:(nrow(tData))) {
      if(sum(!is.na(tData[pep,]))>1 & sum(!is.na(trefData[pep,]))>1) {
        tptvalues<-append(tptvalues,t.test(as.vector(tData[pep,]),as.vector(trefData[pep,]))$p.value)
      } else {
        tptvalues <- append(tptvalues,NA)
      }
    }
    names(tptvalues)<-rownames(tData)
    ptvalues <- cbind(ptvalues,tptvalues)
    
    ## rank products
    # calculate 10 random pairing combinations and then take average of p-values
    NumPairs <- 10
    tpRPvalues<-matrix(NA,ncol=NumPairs,nrow=nrow(Data),dimnames=list(rows = rownames(Data), cols=1:NumPairs))
    for (p in 1:NumPairs) {
      tRPMAData <- tData[,sample(1:NumReps)] - trefData[,sample(1:NumReps)]
      #Up
      RPMAUp_pvalues <- RPStats(tRPMAData,NumReps)
      #Down
      RPMADown_pvalues <- RPStats(-tRPMAData,NumReps)
      ttt <- rowMins(cbind(RPMAUp_pvalues,RPMADown_pvalues),na.rm=T)*2
      ttt[ttt>1] <- 1
      tpRPvalues[names(RPMAUp_pvalues),p] <- ttt
    }
    pRPvalues[,vs-1] <- rowMeans(tpRPvalues,na.rm=T)
    
    ## Permutation tests: add columns from randomized full set to reach min. NumPermCols replicates
    # randomizing also sign to avoid tendencies to one or the other side
    if (ncol(tData)*2<NumPermCols) {
      AddDat <- matrix(sample(as.vector(unlist(tData)),(NumPermCols*0.5-ncol(tData))*nrow(tData),replace=T),nrow=nrow(tData))
      PermData <- cbind(tData, AddDat)
      AddDat <- matrix(sample(as.vector(unlist(trefData)),(NumPermCols*0.5-ncol(trefData))*nrow(trefData),replace=T),nrow=nrow(trefData))
      PermFullData <- cbind(PermData,trefData, AddDat)
    } else {
      PermFullData <- cbind(tData,trefData)
    }
    RealStats <- StatsForPermutTest(cbind(trefData,tData),Paired=F)
    # print(head(PermFullData))
    PermutOut <- parallel::mclapply(1:NTests,function (x) {indat <- apply(PermFullData,1,function(y) sample(y,NumReps*2)*sample(c(1,-1),NumReps*2,replace=T))
    StatsForPermutTest(t(indat),F)
    },mc.cores=NumThreads)
    PermutOut <- matrix(unlist(PermutOut),nrow=nrow(tData))
    PermutOut[!is.finite(PermutOut)] <- NA
    RealStats[!is.finite(RealStats)] <- NA
    pPermutvalues[,vs-1] <- apply(cbind(RealStats,PermutOut), 1 , function(x) ifelse(is.na(x[1]) | sum(!is.na(x)) == 0,NA,(1+sum(x[1] < x[-1],na.rm=T))/(sum(!is.na(x)))))
    # print(PermMAData[1,])
    # print(pPermutvalues)
  }
  lratios <- NULL
  qRPvalues <- qtvalues <- qPermutvalues <- matrix(NA,nrow=nrow(Data),ncol=NumCond-1,dimnames=list(rows=rownames(Data), cols=1:(NumCond-1)))
  for (i in 1:(NumCond-1)) {
    tqs <- qvalue(na.omit(ptvalues[,i]))$qvalues
    qtvalues[names(tqs),i] <- tqs
    print(range(pPermutvalues[,i]))
    tqs <- qvalue(na.omit(pPermutvalues[,i]))$qvalues
    qPermutvalues[names(tqs),i] <- tqs
    print(range(na.omit(pRPvalues[,i])))
    tqs <- qvalue(na.omit(pRPvalues[,i]),lambda=seq(0.05,max(na.omit(pRPvalues[,i]))-0.05,0.05))$qvalues
    qRPvalues[names(tqs),i] <- tqs
    lratios <- cbind(lratios, rowMeans(Data[,Reps==i+1],na.rm=T)-rowMeans(Data[,Reps==1],na.rm=T))
    
  }
  
  return(list(lratios=lratios,ptvalues=ptvalues, plvalues=plvalues, pRPvalues=pRPvalues, pPermutvalues=pPermutvalues,
              qtvalues=qtvalues, qlvalues=qlvalues, qRPvalues=qRPvalues,qPermutvalues=qPermutvalues,Sds=sqrt(lm.bayes$s2.post)))
}

MissingStats <- function(Data, NumCond, NumReps) {
  Reps <- rep(1:NumCond,NumReps)
  pNAvalues<-matrix(NA,ncol=NumCond-1,nrow=nrow(Data),dimnames=list(rows = rownames(Data), cols=1:(NumCond-1)))
  qNAvalues<-matrix(NA,ncol=NumCond-1,nrow=nrow(Data),dimnames=list(rows = rownames(Data), cols=1:(NumCond-1)))
  for (vs in 2:NumCond) {
    tData<-Data[,Reps==vs]
    trefData <- Data[,Reps==1]
    tCompDat <- cbind(tData,trefData)
    qs <- quantile(tCompDat, probs=seq(0,1,0.01),na.rm=T)
    pvals <- statis <- matrix(NA,nrow(tCompDat),ncol=length(qs))
    for (q in qs) {
      tCompDat[tCompDat<q] <- NA
      # print(head(tCompDat))
      NAPDistr <- MissValPDistr(NumReps,sum(is.na(tCompDat))/(nrow(tCompDat)*2*NumReps))
      statis[,which(q==qs)] <- (rowSums(!is.na(tCompDat[,1:NumReps])) - rowSums(!is.na(tCompDat[,(NumReps+1):(2*NumReps)])))
      pvals[,which(q==qs)] <- NAPDistr[abs(statis[,which(q==qs)])+1]
      
    }
    pNAvalues[,vs-1] <- rowMins(pvals)*(NumReps+1)
    # qNAvalues[,vs-1] <- qvalue(pNAvalues[,vs-1],lambda=seq(0.1,max(pNAvalues[,vs-1]),length=100))$qvalue
    qNAvalues[,vs-1] <- p.adjust(pNAvalues[,vs-1], method="BH")
    
  }
  head(pNAvalues)
  pNAvalues[pNAvalues>1] <- 1
  
  return(list(pNAvalues=pNAvalues, qNAvalues=qNAvalues))
  
}

# Function to determine "optimal" fold-change and q-value thresholds
# the idea is to maximize the percental output of features commonly found for limma, rank products, permutation and NA tests
FindFCandQlim <- function(Qvalue, LogRatios, NumTests) {
  
  BestComb <- c(0,0)
  BestRegs <- 0
  NumCond <- ncol(LogRatios)+1
  
  Qvalue[is.na(Qvalue)] <- 1
  
  smallestq <- signif(min(Qvalue,na.rm=T))
  qrange <- c(0.1,0.2,0.5)*10^(rep(-10:0,each=3))
  qrange <- qrange[which.min(abs(smallestq-qrange)):(length(qrange)-2)]
  
  # Run over different FC thresholds 
  fcRange <- seq(0,max(abs(range(LogRatios,na.rm=T))),length=100)
  BestVals <- parallel::mclapply(fcRange, function(fc) { 
    for (t in 1:(NumTests-2)) {
      tvals <- Qvalue[,(NumCond-1)*t+1:(NumCond-1)]
      tvals[LogRatios < fc & LogRatios > -fc] <- 1
      Qvalue[,(NumCond-1)*t+1:(NumCond-1)] <- tvals
    }
    # Run over range of q-values
    for (qlim in qrange) {
      alldistr <- vector("numeric",NumCond-1)
      for (t in 1:(NumCond-1)) {
        distr <- table(rowSums(Qvalue[,seq(t,ncol(Qvalue),NumCond-1)] < qlim, na.rm=T))
        allregs <- sum(distr[2:length(distr)],na.rm=T)
        # print(distr)
        if (length(distr) > 1 & !is.na(distr["4"]))
          alldistr[t] <- distr["4"]/allregs
      }
      # print(alldistr)
      if (mean(alldistr) > BestRegs) {
        BestRegs <- mean(alldistr)
        BestComb <- c(fc,qlim)
        # print(BestRegs)
        
      }
    }
    c(BestRegs,BestComb)
  },mc.cores=NumThreads)
  
  BestRegs <- 0
  for (i in 1:length(BestVals)) {
    if (BestRegs < BestVals[[i]][1]){
      BestComb <- BestVals[[i]][2:3]
      BestRegs <- BestVals[[i]][1]
    }
    # print(BestVals)
  }
  return(BestComb)
  
}


# calculated common q-value over different tests. Hommel methods gives 
# upper bound for p-values coming from independent or positively dependent tests
UnifyQvals <- function(Qvalue, NumCond, NumTests) {
  UnifiedQvalue <- matrix(NA,ncol=NumCond-1,nrow=nrow(Qvalue))
  for (i in 1:(NumCond-1)) {
    # print(seq(i,ncol(Qvalue)-NumCond+1,NumCond-1))
    UnifiedQvalue[,i] <- colMins(apply(Qvalue[,seq(i,ncol(Qvalue)-NumCond+1,NumCond-1)], 1, p.adjust, "hommel"),na.rm=T)
  }
  UnifiedQvalue
  
  
}

