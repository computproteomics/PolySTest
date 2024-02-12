# Function to set mfrow to a maximum of max_col columns
set_mfrow <- function(num_total, max_col) {
  if (num_total <= max_col) {
    par(mfrow = c(1, num_total))
  } else {
    par(mfrow = c(ceiling(num_total / max_col), max_col))
  }
}

#' plotPvalueDistr
#'
#' Function to plot p-value distributions for all tests and comparisons
plotPvalueDistr <- function(fulldata, testNames, compNames,
                            testCols=c("#33AAAA","#33AA33","#AA3333","#AA33AA",
                                       "#AAAA33","#3333AA"), ...) {
  print("Plot p-values")
  rdat <- rowData(fulldata)
  PValue <- as.matrix(rdat[,grep("^p_value", colnames(rdat)),drop=F])
  NumTests <- length(testNames)
  NumComps <- length(compNames)
  set_mfrow((NumTests)*NumComps, NumTests)
  if (ncol(PValue) > 0) {
    for (i in seq_len(NumComps)) {
      for (j in 1:NumTests) {
        hist(PValue[,(NumComps)*(j-1)+i],100, main=testNames[j],
             sub=compNames[i],col=testCols[j+1],
             xlab="p-value",border=NA, ...)
      }
    }
  }
}

#' plotVolcano
#'
#' Volcano plot for all tests and all comparisons
plotVolcano <- function(fulldata, testNames, compNames, sel_prots, qlim=0.05, fclim=c(0,0),
                        testCols=c("#33AAAA","#33AA33","#AA3333","#AA33AA",
                                   "#AAAA33","#3333AA")) {
  print("Plot volcano plots")
  rdat <- rowData(fulldata)
  Qvalue <- as.matrix(rdat[,grep("^FDR", colnames(rdat)),drop=F])
  if (ncol(Qvalue) > 0) {
    LogRatios <- as.matrix(rdat[,grep("^log_ratios", colnames(rdat))])
    NumTests <- length(testNames)
    NumComps <- length(compNames)
    set_mfrow(NumTests*NumComps, NumTests)
    for (i in 1:(NumComps)) {
      for (j in 1:NumTests) {
        plot(LogRatios[,i],-log10(Qvalue[,(NumComps)*(j-1)+i]),
             main=testNames[j],sub=compNames[i],
             xlab="log fold-change",ylab="-log10(q)",
             cex=colSelected(0.5,nrow(Qvalue),sel_prots,1),
             col=colSelected(adjustcolor(testCols[j],alpha.f=0.3),nrow(Qvalue),
                             sel_prots,"#FF9933"),pch=16,
             ylim=-log10(c(1,min(Qvalue,na.rm=T))))

        abline(h=-log10(qlim),col="#AA3333",lwd=2)
        abline(v=fclim,col="#AA3333",lwd=2)
      }
    }
  }
}


plotExpression <- function(fulldata, sel_prots, testNames, compNames, profiles_scale=TRUE, qlim=0.05, fclim=c(0,0)) {
  NumComps <- length(compNames)
  print(compNames)
  NumTests <- length(testNames)
  rdat <- rowData(fulldata)
  LogRatios <- as.matrix(rdat[,grep("^log_ratios", colnames(rdat)),drop=F])
  Qvalue <- as.matrix(rdat[,grep("^FDR", colnames(rdat)),drop=F])
  rownames(Qvalue) <- rownames(LogRatios) <-  rownames(rdat)
  if (ncol(Qvalue) > 0 & length(sel_prots) > 1) {
    FCRegs <- filterFC(fulldata, NumTests, NumComps, fclim)
    NumCond <- metadata(fulldata)$NumCond
    NumReps <- metadata(fulldata)$NumReps

    # CI plots of max 30 features
    print(sel_prots)
    print(dim(Qvalue))
    SubSetQval <- Qvalue[sel_prots,,drop=F]
    SubSetLR <- LogRatios[sel_prots,,drop=F]
    SubSetLR <- SubSetLR[order(rowMins(SubSetQval[,1:(NumComps),drop=F],na.rm=T)),,drop=F]

    SubSet <- SubSetLR[1:min(nrow(SubSetLR),30),,drop=F]
    indices <- rownames(SubSet)
    tdat <- as.matrix(dat[rownames(SubSet), (rep(1:NumReps,NumCond)-1)*NumCond+rep(1:NumCond,each=NumReps),drop=F])
    rownames(tdat) <- strtrim(rownames(tdat), 20)
    MeanSet <- SDSet <-  matrix(NA,nrow=nrow(tdat),ncol=NumCond,dimnames =
                                  list(x=rownames(tdat),y=paste("Condition",1:NumCond)))
    for (c in 1:NumCond) {
      MeanSet[,c] <- rowMeans(tdat[,1:NumReps + (c-1)*NumReps,drop=F],na.rm=T)

      SDSet[,c] <- rowSds(tdat[,1:NumReps + (c-1)*NumReps,drop=F],na.rm=T)
    }
    if (profiles_scale) {
      MeanSet <- MeanSet - rowMeans(MeanSet,na.rm=T)
    }

    layout(t(c(1,1,2,2,3,3)))
    plot(0,0,type="n",bty="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA)
    legend("topright",col=rainbow(nrow(SubSet),alpha = 0.8,s=0.7),legend=strtrim(rownames(SubSet),20),lwd=3,
           title="Features")
    plotCI(1:(NumCond)+runif(1,-0.1,0.1),MeanSet[1,],pch=16,
           xlab="Conditions",xlim=c(0.7,(NumCond+1)-0.7),
           ylab="expression values",col=rainbow(nrow(MeanSet),alpha = 0.8,s=0.7)[1],
           uiw=SDSet[1,],type="b",barcol="#000000AA",
           ylim=range(MeanSet,na.rm=T),xaxt="none",lwd=1.5)
    title(main="Feature expression over conditions")
    axis(1,at=1:(NumCond),labels = colnames(MeanSet))
    # abline(h=fclim)
    if (nrow(MeanSet)>1) {
      for (i in 2:nrow(MeanSet)){
        plotCI(1:(NumCond)+runif(1,-0.1,0.1),MeanSet[i,],
               add = T,pch=16,col=rainbow(nrow(MeanSet))[i],
               uiw=SDSet[i,],type="b",barcol="#000000AA",lwd=1.5)


      }
    }

    if (length(SubSet)> 0) {
      par(mar=rep(0,4))
      circos.clear()
      circos.par(cell.padding=c(0,0,0,0),canvas.xlim=c(-1.5,1.5),canvas.ylim=c(-1.5,1.5),
                 track.margin=c(0,0.02),start.degree=90,gap.degree=4)
      circos.initialize(1:(NumComps),xlim=c(0,1))
      for (t in 1:NumTests) {
        # print(tsign)
        nfeat <- min(nrow(SubSet),30)
        cols <- rainbow(nfeat,alpha = 0.8,s=0.7)
        tsign <- FCRegs[indices,(t-1)*(NumComps)+(1:(NumComps)),drop=F]<qlim
        circos.trackPlotRegion(ylim=c(-3,2),track.height=1/12, bg.border="#777777",
                               panel.fun = function(x, y) {
                                 name = get.cell.meta.data("sector.index")
                                 i = get.cell.meta.data("sector.numeric.index")
                                 xlim = get.cell.meta.data("xlim")
                                 ylim = get.cell.meta.data("ylim")
                                 xdiff <- (xlim[2]-xlim[1])/nfeat
                                 if(t == 1) {
                                   circos.text(mean(xlim), max(ylim)+30, compNames[i], facing = "inside",
                                               niceFacing = TRUE,cex = 1,font=2)
                                   circos.axis("top", labels = strtrim(rownames(SubSetLR),20),
                                               major.at=seq(1/(nfeat*2),1-1/(nfeat*2),length=nfeat),minor.ticks=0,
                                               labels.cex = 0.8,labels.facing = "reverse.clockwise",
                                   )
                                 }
                                 for (j in which(tsign[,i])) {
                                   circos.rect(xleft=xlim[1]+(j-1)*xdiff, ybottom=ylim[1],
                                               xright=xlim[2]-(nfeat-j)*xdiff, ytop=ylim[2],
                                               col = cols[j], border=NA)
                                 }
                               })
      }
      fccols <- redblue(1001)
      # print(SubSet)
      circos.trackPlotRegion(ylim=c(-3,2),track.height=1/4, bg.border=NA, panel.fun = function(x, y) {
        name = get.cell.meta.data("sector.index")
        i = get.cell.meta.data("sector.numeric.index")
        xlim = get.cell.meta.data("xlim")
        ylim = get.cell.meta.data("ylim")
        #circos.text(x=mean(xlim), y=1.7,
        #            labels=name, facing = dd, cex=0.6,  adj = aa),
        xdiff <- (xlim[2]-xlim[1])/nfeat
        for (j in 1:nfeat) {
          # print((SubSetLR[j,i]/max(LogRatios,na.rm=T))*500+500)
          circos.rect(xleft=xlim[1]+(j-1)*xdiff, ybottom=ylim[1],
                      xright=xlim[2]-(nfeat-j)*xdiff, ytop=ylim[2],
                      col = fccols[(SubSetLR[j,i]/max(LogRatios,na.rm=T))*500+500], border=0)
        }})
      text(0,0,"Log\nratios",cex=0.7)
      # label the different tracks
      mtext(paste("Successful statistical tests\nfor threshold given above.\nFrom outer to inner circles",
                  paste("Track ",1:NumTests,": ",testNames,sep="",collapse="\n"),sep="\n"),
            side=1,outer=T,adj=1,line=-1,cex=0.6)
    }
  }
}

plotUpset <- function (fulldata, NumTests, NumComps, qlim=0.05, fclim=c(0,0)) {

  FCRegs <- filterFC(fulldata,  NumTests, NumComps, fclim)

  if (!is.null(FCRegs)) {

    WhereRegs <- FCRegs[,rep(0:(NumTests-2), NumComps)*(NumComps)+rep(1:(NumComps),each=NumTests-1),drop=F]<qlim
    WhereRegs[WhereRegs] <- 1
    deleted_cols <- which(colSums(WhereRegs,na.rm=T)==0)
    # print(deleted_cols)

    tcolnames <- paste("A",rep(1:(NumComps),each=NumTests-1))
    if (length(deleted_cols) > 0) {
      tcolnames <- tcolnames[-deleted_cols]
      WhereRegs <- WhereRegs[,-deleted_cols,drop=F]
    }
    tcols <- rep(rainbow(NumComps),each=1)
    names(tcols) = rep(paste("A",1:(NumComps)),1)


    if(length(WhereRegs)>0) {
      upset(as.data.frame(WhereRegs),nsets=ncol(WhereRegs),mainbar.y.label = "Significant features",
            nintersects = NA,keep.order=F,sets=colnames(WhereRegs),text.scale=1.5, mb.ratio = c(0.55, 0.45),
            set.metadata = list(data = data.frame(set=colnames(WhereRegs),cols=tcolnames,crab=1:ncol(WhereRegs)),
                                plots = list(list(type = "matrix_rows",column = "cols", colors=tcols,alpha=0.5))))
    } else {
      NULL
    }
  }
}


plotRegDistr <- function( fulldata, testNames, NumComps, qlim=0.05, fclim=c(0,0),
                          TestCols=c("#33AAAA","#33AA33","#AA3333","#AA33AA",
                                     "#AAAA33","#3333AA")) {
  NumTests <- length(testNames)
  FCRegs <- filterFC(fulldata, NumTests, NumComps, fclim)

  if (!is.null(FCRegs)) {

    set_mfrow(NumComps, 5
              )

    rdat <- rowData(fulldata)
    Qvalue <- as.matrix(rdat[,grep("^FDR", colnames(rdat)),drop=F])
    tmpX <- 10^seq(log10(min(Qvalue,na.rm=T)),0.1,0.01)

    for (i in 1:(NumComps)) {
      plot(tmpX,rowSums(sapply((FCRegs[,i]),"<",tmpX),na.rm=T), main=paste("Comparison",i),xlab="FDR threshold",
           ylab="Number significant",type="l",col=TestCols[1],ylim=c(1,nrow(Qvalue)),log="xy",lwd=2)
      if (i==1)
        legend("topleft",legend = testNames,col=TestCols,lwd=2)
      lines(tmpX,rowSums(sapply((FCRegs[,(NumComps)*1+i]),"<",tmpX),na.rm=T),col=TestCols[2],lwd=2)
      lines(tmpX,rowSums(sapply((FCRegs[,(NumComps)*2+i]),"<",tmpX),na.rm=T),col=TestCols[3],lwd=2)
      lines(tmpX,rowSums(sapply((FCRegs[,(NumComps)*3+i]),"<",tmpX),na.rm=T),col=TestCols[4],lwd=2)
      lines(tmpX,rowSums(sapply((FCRegs[,(NumComps)*4+i]),"<",tmpX),na.rm=T),col=TestCols[5],lwd=2)
      lines(tmpX,rowSums(sapply((FCRegs[,(NumComps)*5+i]),"<",tmpX),na.rm=T),col=TestCols[6],lwd=2)
      abline(v=qlim,col=2)
    }
  }
}


plotHeatmaply <- function( fulldata, sel_prots, NumComps, heatmap_scale="none", file=NULL) {
  p <- plotly_empty()

  rdat <- rowData(fulldata)
  LogRatios <- as.matrix(rdat[,grep("^log_ratios", colnames(rdat)),drop=F])
  Qvalue <- as.matrix(rdat[,grep("^FDR", colnames(rdat)),drop=F])
  rownames(LogRatios) <- rownames(Qvalue) <- rownames(rdat)
  if (ncol(Qvalue) > 0 & length(sel_prots) > 0) {

    NumCond <- metadata(fulldata)$NumCond
    NumReps <- metadata(fulldata)$NumReps

    # CI plots of max 30 features
    SubSetQval <- Qvalue[sel_prots,,drop=F]
    SubSetLR <- LogRatios[sel_prots,,drop=F]
    SubSetLR <- SubSetLR[order(rowMins(SubSetQval[,1:(NumComps),drop=F],na.rm=T)),,drop=F]

    if (!is.null(SubSetLR)) {
      if (length(SubSetLR)> 0 & nrow(SubSetLR)>1) {
        print("running heatmap")
        withProgress(message="Creating heatmap ...", min=0,max=1, {
          setProgress(0.5)
          tdat <- dat[rownames(SubSetLR), (rep(1:NumReps,NumCond)-1)*NumCond+rep(1:NumCond,each=NumReps),drop=F]
          rownames(tdat) <- strtrim(rownames(tdat), 30)

          # remove data rows with more than 45% missing values
          to_remove <- which(rowSums(is.na(tdat)) > ncol(tdat)*0.45)
          tqvals <- Qvalue[rownames(SubSetLR),1:(NumComps),drop=F]
          if (length(to_remove)>0) {
            tqvals <- tqvals[-to_remove,,drop=F]
            tdat <- tdat[-to_remove,,drop=F]
          }
          # setting colors of p-values
          pcols <- rev(c(0.001, 0.01, 0.05, 1))
          ttt <- tqvals
          for (c in pcols) {
            ttt[tqvals <= c] <- c
          }
          tqvals <- data.frame(ttt)
          for (c in 1:ncol(tqvals))
            tqvals[,c] <- paste("<",as.character(tqvals[,c],pcols),sep="")
          scaling <- "none"
          tqvals <- tqvals[order(rownames(tdat)),]
          if(heatmap_scale) scaling <- "row"
          p <- heatmaply(tdat[order(rownames(tdat)),,drop=F],Colv=F,scale =scaling,trace="none",cexRow=0.7,plot_method="plotly",
                         RowSideColors = tqvals, row_side_palette = grey.colors, file=file)
        })
      }

    }
  }
  p
}
