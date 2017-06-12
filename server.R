# Visualization of expression patterns + q-values comparison (of the marked ones) + usetR package for comparison of what is regulated (select directly them to display)
# BUtton to deselect all features. Click on figures to select for expressions?
# In q-q plots: visualize by areas? visual helps? correlations?


# library(shinyBS)
library(DT)
library(circlize)
library(UpSetR)
library(matrixStats)
library(qvalue)
library(gplots)
source("HelperFuncs.R")
shinyServer(function(input, output,clientData,session) {
  dat <- NULL
  addInfo <- NULL # additional info to be re-added after analysis
  heightSize = function() {
    (input$NumCond-1)*200
  }
  heightSize2 = function() {
    (input$NumCond-1)*400
  }
  colSelected = function(col, num, sel, col2) {
    ttt <- rep(col,num)
    ttt[sel] <- col2
    return(ttt)
  }
  
  output$messages <- renderText("")
  output$thead <- renderText("Multiple statistical tests for the detection for differentially regulated features")
  output$description <- renderText("The methods are designed to carry out statistical tests on data with few replicates and 
                                   high amounts of missing values. ")
  output$description2 <- renderText("We recommend this method for datasets with a minimum of 1000 features (rows of the table)
                                   and at least 3 replicates. The method is based on ratios between
features within the replicate, i.e. the tests are carried out on paired tests.")
  output$description3 <- renderText("For details, see <i>in preparation</i>
                                  ")
  output$description4 <- renderText("Data input: Data table (csv file) with separate columns for feature information and quantification values. If row names is set,
                                    then the first column of the csv file provides unique features names (no duplicates).
                                    Columns with quantifications follow after all other columns and start with columns number <i>First column for quantification</i>.The quantification values should be 
                                    log-transformed intensity/abundance values (not ratios) which have been normalized to be comparable. The 
                                   order of the columns is required to be A1, A2, A3, ..., B1, B2, B3, ..., where
                                   1,2,3 ... are the conditions and A,B,... denote replicates.
The tests check for differentially regulated features
                                    versus the \"reference\" condition. For each comparison (log-ratio), we plot the histograms of the uncorrected p-values
                                    and volcano plots of the false discovery rates 
                                    (corrected for multiple testing according to Storey JD. All points above the horizontal line have a q-value below 0.05. A direct approach to false discovery rates. Journal of the Royal Statistical Society. 2002;64:479–498).")
  #   output$NumReps <- renderText("Number of replicates")
  #   output$NumCond <- renderText("Number of conditions")
  #   output$refCond <- renderText("Reference condition")
  
  
  output$downloadData <- downloadHandler(
    filename = function() {
      validate(
        need(F,"No data")
      )      
      paste("Results", Sys.Date(), ".csv", sep="");
    },
    content = function(file) {
      write.csv(FullReg, file)
    }
  )
  
  observeEvent(input$example, {
    dat <<- "EXAMPLE"
    # TODO update first column slider
    updateSliderInput(session,"NumCond",value=4)
    updateSliderInput(session,"NumReps",value=3)
    updateSliderInput(session,"refCond",value=1)
  })
  
  output$plotpval <- renderPlot({
    input$button
    # input$example
    dev.control(displaylist="enable")
    updateNumericInput(session,"refCond",max=input$NumCond,value=input$refCond)
    output$messages <- renderText("")
    if(!is.null(dat)) {
      if (dat != "EXAMPLE")
        dat <- input$in_file
    } else {
      dat <- input$in_file
    }
    NumReps <- input$NumReps
    NumCond <- input$NumCond
    refCond <- input$refCond
    isPaired <- input$is_paired
    # print(dat)
    
    if (!is.null(dat)) {
      if (dat == "EXAMPLE")  {
        dat <- read.csv("LiverAllProteins.csv",row.names=1)
      } else {
        delim <- input$delimiter
        if (delim == "tab")
          delim <- "\t"
        if (input$row.names){
          dat <- read.csv(input$in_file$datapath,row.names=1,header=input$is_header,sep=delim,dec=input$digits)
          updateSliderInput(session,"QuantCol",max=ncol(dat))
          
          # delete row with empty name
          dat <- dat[!rownames(dat)=="",]
          if (input$ColQuant > 2) {
            addInfo <- dat[,1:(input$ColQuant-2),drop=F]
            for (c in 1:ncol(addInfo))
              addInfo[,c] <- as.character(addInfo[,c])
            print(head(addInfo))
            dat <- dat[,-(1:(input$ColQuant-2))]
          }
        } else {
          dat <- read.csv(input$in_file$datapath,header=input$is_header,sep=delim,dec=input$digits)
          updateSliderInput(session,"QuantCol",max=ncol(dat))
          
          if (input$ColQuant > 1) {
            addInfo <- dat[,1:(input$ColQuant-1),drop=F]
            for (c in 1:ncol(addInfo))
              addInfo[,c] <- as.character(addInfo[,c])
            
            dat <- dat[,-(1:(input$ColQuant-1))]
          }
        }
      }
      # print(head(dat))
      updateSliderInput(session,"NumCond",max=ncol(dat))
      updateSliderInput(session,"NumReps",max=ncol(dat))
      updateSliderInput(session,"refCond",max=input$NumCond)
      
      output$input_stats <- renderText(paste(ifelse(mode(as.matrix(dat))!="numeric","<b>Wrong file format /setup</b></br>",""),
                                             ifelse(ncol(dat) != NumReps*NumCond,"<b>Column number doesn't fit with number of replicates and conditions!</b><br/>",""),
                                             "Number of features: ",nrow(dat),
                                             "<br/>Percentage of missing values:",
                                             round(sum(is.na(dat))/nrow(dat)/ncol(dat)*100,digits = 2),"<br/>",
                                             paste("<br/>Comparison",1:(NumCond-1),": Condition",
                                                   (1:NumCond)[-(input$refCond)],"versus",input$refCond,collapse=""),"<br/>"))
      isolate({
        if (input$button == 0)
          return()
        if (ncol(dat) == NumReps*NumCond) {
          withProgress(message="Calculating ...", min=0,max=1, {
            
            #           output$messages <- renderText("running")
            incProgress(0.1, detail = paste("Running statistical tests"))
            dat [!is.finite(as.matrix(dat))] <- NA
            
            Against<-seq(refCond,NumCond*NumReps,NumCond)
            RR<-1:(NumReps*NumCond)
            RR<-RR[-Against]
            RR<-rbind(RR,rep(Against,each=NumCond-1))
            MAData<-NULL
            # Rearranged dat (reference condition comes first)
            UData <- NULL
            for (i in 1:NumReps) {
              UData<-cbind(UData,dat[,(NumCond)*(i-1)+refCond])
              UData <- cbind(UData,dat[,(NumCond)*(i-1)+(1:NumCond)[-refCond]])
            }
            rownames(UData) <- rownames(dat)
            if (isPaired) {
              for (i in 1:ncol(RR)) {
                MAData<-cbind(MAData,dat[,RR[1,i]]-dat[,RR[2,i]])
              }
              rownames(MAData)<-rownames(dat)
              qvalues <- Paired(MAData, NumCond-1, NumReps)
            } else {
              qvalues <- Unpaired(UData, NumCond, NumReps)
            }
            # get q-values for missing values
            MissingStats <- MissingStats(UData, NumCond, NumReps)
            
            
            incProgress(0.5, detail = paste("Preparing data"))
            
            LogRatios <- qvalues$lratios
            Pvalue <- cbind(qvalues$ptvalues, qvalues$plvalues, qvalues$pRPvalues, qvalues$pPermutvalues, MissingStats$pNAvalues)
            Qvalue <- cbind(qvalues$qtvalues, qvalues$qlvalues, qvalues$qRPvalues, qvalues$qPermutvalues, MissingStats$qNAvalues)
            # WhereReg <- cbind(qvalues$qtvalues<qlim, qvalues$qlvalues<qlim, qvalues$qRPvalues<qlim, qvalues$qPermutvalues<qlim, MissingStats$qNAvalues<qlim)
            
            compNames <- paste("C",RR[1,1:(NumCond-1)]," vs C",RR[2,1:(NumCond-1)],sep="")
            testNames <- c("t-test","limma","rank products","Permutation test","NA test")
            colnames(LogRatios) <- paste("log-ratios",compNames)
            colnames(Pvalue) <- paste("p-values",rep(testNames,each=NumCond-1),rep(compNames,length(testNames)))
            colnames(Qvalue) <- paste("q-values",rep(testNames,each=NumCond-1),rep(compNames,length(testNames)))
            # colnames(WhereReg) <- paste("Differentially regulated",rep(testNames,each=NumCond-1),rep(compNames,length(testNames)))
            
            # print(cor(log10(Pvalue),use="na.or.complete"))
            if (!is.null(addInfo))
              FullReg <- cbind(addInfo[rownames(LogRatios),],LogRatios, Qvalue)#, WhereReg)
            else 
              FullReg <- cbind(LogRatios, Qvalue)#, WhereReg)
            
            output$stat_table <- DT::renderDataTable(FullReg,server=T,extensions=c("Buttons","FixedColumns"),
                                                     options = list(scrollX = TRUE,dom = 'Blfrtip',#fixedColumns=T,
                                                                    lengthMenu = c(5, 10, 50,100),
                                                                    buttons = list('colvis', 'copy', 'print' 
                                                                    )))
            
            proxy = dataTableProxy('stat_table')
            
            observeEvent(input$resetSelection, {
              proxy %>% selectRows(NULL)
            })
            
            observeEvent(input$plotregdistr_click, {
              print(input$plotregdistr_click$x)
              print(input$plotregdistr_click$y)
            })
            
            
            
            
            incProgress(0.7, detail = paste("Plotting first results"))
            
            
            ## Plotting DRFs vs q-value thresholds + volcano plots + p-val vs q-vals
            par(mfrow=c(NumCond-1,length(testNames)))
            for (i in 1:(NumCond-1)) {
              hist(Pvalue[,i],100, main="t-test",col="#333366",xlab="p-value")
              hist(Pvalue[,(NumCond-1)+i],100, main=paste(compNames[i],"\nlimma",sep="\n"),col="#336633",xlab="p-value")
              hist(Pvalue[,(NumCond-1)*2+i],100, main="rank products",col="#663333",xlab="p-value")
              hist(Pvalue[,(NumCond-1)*3+i],100, main="Permutation test",col="#663366",xlab="p-value")
              hist(Pvalue[,(NumCond-1)*4+i],100, main="NA test",col="#666633",xlab="p-value")
            }
            par(mfrow=c(1,1))
            # volcano plots
            output$plotvolc <- renderPlot({
              qlim <- input$qval
              input$button
              input$stat_table
              par(mfrow=c(NumCond-1,5))
              for (i in 1:(NumCond-1)) {
                plot(LogRatios[,i],-log10(Qvalue[,i]), main="t-test",xlab="log fold-change",ylab="-log10(q)",
                     cex=colSelected(1,nrow(Qvalue),input$stat_table_rows_selected,2),
                     col=colSelected("#33336655",nrow(Qvalue),input$stat_table_rows_selected,"#FF9933"),pch=16,ylim=-log10(c(1,min(Qvalue,na.rm=T))))
                abline(h=-log10(qlim),col="#AA3333",lwd=2)
                plot(LogRatios[,i],-log10(Qvalue[,(NumCond-1)+i]), main=paste(compNames[i],"limma",sep="\n"),xlab="log fold-change",ylab="-log10(q)",
                     cex=colSelected(1,nrow(Qvalue),input$stat_table_rows_selected,2),
                     col=colSelected("#33663355",nrow(Qvalue),input$stat_table_rows_selected,"#FF9933"),pch=16,ylim=-log10(c(1,min(Qvalue,na.rm=T))))
                abline(h=-log10(qlim),col="#AA3333",lwd=2)
                plot(LogRatios[,i],-log10(Qvalue[,(NumCond-1)*2+i]), main="rank products",xlab="log fold-change",ylab="-log10(q)",
                     cex=colSelected(1,nrow(Qvalue),input$stat_table_rows_selected,2),
                     col=colSelected("#66333355",nrow(Qvalue),input$stat_table_rows_selected,"#FF9933"),pch=16,ylim=-log10(c(1,min(Qvalue,na.rm=T))))
                abline(h=-log10(qlim),col="#AA3333",lwd=2)
                plot(LogRatios[,i],-log10(Qvalue[,(NumCond-1)*3+i]), main="Permutation tests",xlab="log fold-change",ylab="-log10(q)",
                     cex=colSelected(1,nrow(Qvalue),input$stat_table_rows_selected,2),
                     col=colSelected("#66336655",nrow(Qvalue),input$stat_table_rows_selected,"#FF9933"),pch=16,ylim=-log10(c(1,min(Qvalue,na.rm=T))))
                abline(h=-log10(qlim),col="#AA3333",lwd=2)
                plot(LogRatios[,i],-log10(Qvalue[,(NumCond-1)*4+i]), main="NA tests",xlab="log fold-change",ylab="-log10(q)",
                     cex=colSelected(1,nrow(Qvalue),input$stat_table_rows_selected,2),
                     col=colSelected("#66663355",nrow(Qvalue),input$stat_table_rows_selected,"#FF9933"),pch=16,ylim=-log10(c(1,min(Qvalue,na.rm=T))))
                abline(h=-log10(qlim),col="#AA3333",lwd=2)
              }
              par(mfrow=c(1,1))
              
            },height=heightSize2)
            
            output$plotexpression <- renderPlot({
              input$button
              input$stat_table
              qlim <- input$qval
              indices <- input$stat_table_rows_selected
              SubSet <- FullReg[input$stat_table_rows_selected,,drop=F]
              if (length(SubSet)> 0) {
                par(mfrow=c(2,2))
                plot(0,0,type="n",bty="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA)
                # print(colnames(SubSet))
                legend("topright",col=rainbow(nrow(SubSet),alpha = 0.8,s=0.7),legend=rownames(SubSet),lwd=3)
                plotCI(1:(NumCond-1)+runif(1,-0.1,0.1),LogRatios[as.vector(indices)[1],1:(NumCond-1),drop=F],pch=16,
                       xlab="Conditions",xlim=c(0.7,(NumCond)-0.7),
                       ylab="log-ratios",col=rainbow(nrow(SubSet),alpha = 0.8,s=0.7)[1],
                       uiw=qvalues$Sds[input$stat_table_rows_selected][1],type="b",barcol="#000000FF",
                       ylim=range(LogRatios,na.rm=T),xaxt="none",lwd=1.5)
                axis(1,at=1:(NumCond-1),labels = compNames)
                if (nrow(SubSet)>1) {
                  for (i in 2:min(nrow(SubSet),20)){
                    plotCI(1:(NumCond-1)+runif(1,-0.1,0.1),LogRatios[as.vector(indices[i]),1:(NumCond-1),drop=F],
                           add = T,pch=16,col=rainbow(nrow(SubSet))[i],
                           uiw=qvalues$Sds[input$stat_table_rows_selected][i],type="b",barcol="#000000FF",lwd=1.5) 
                    
                  }
                }
                # Missing values
                # find descent visualization 
                
                # Compare q-values for each protein and method
                SubSet <- Qvalue[input$stat_table_rows_selected,,drop=F]
                # validate(
                #   need(length(SubSet)>0, "Please select features from the data table")
                # )
                
                if (length(SubSet)> 0) {
            
                par(mar=rep(0,4))
                circos.clear()
                circos.par(cell.padding=c(0,0,0,0),canvas.xlim=c(-1.5,1.5),canvas.ylim=c(-1.5,1.5),
                           track.margin=c(0,0.02),start.degree=90,gap.degree=4)
                circos.initialize(1:(NumCond-1),xlim=c(0,1))
                for (t in 1:5) {
                  tsign <- SubSet[,(t-1)*(NumCond-1)+(1:(NumCond-1)),drop=F]<qlim
                  tsign[is.na(tsign)] <- F
                  # print(tsign)
                  nfeat <- min(nrow(SubSet),20)
                  cols <- rainbow(nfeat,alpha = 0.8,s=0.7)

                  circos.trackPlotRegion(ylim=c(-3,2),track.height=1/10, bg.border="#777777", 
                                         panel.fun = function(x, y) {
                                           name = get.cell.meta.data("sector.index")
                                           i = get.cell.meta.data("sector.numeric.index")
                                           xlim = get.cell.meta.data("xlim")
                                           # print(nfeat)
                                           ylim = get.cell.meta.data("ylim")
                                           xdiff <- (xlim[2]-xlim[1])/nfeat
                                           if(t == 1) {
                                             circos.text(mean(xlim), max(ylim)+30, compNames[i], facing = "inside", niceFacing = TRUE,cex = 1,font=2)
                                             circos.axis("top", labels = rownames(SubSet),major.at=seq(1/(nfeat*2),1-1/(nfeat*2),length=nfeat),minor.ticks=0,
                                                         labels.cex = 1,labels.facing = "reverse.clockwise")
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
                SubSet <- LogRatios[input$stat_table_rows_selected,,drop=F]
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
                    print((SubSet[j,i]/max(LogRatios,na.rm=T))*500+500)
                    circos.rect(xleft=xlim[1]+(j-1)*xdiff, ybottom=ylim[1],
                                xright=xlim[2]-(nfeat-j)*xdiff, ytop=ylim[2],
                                col = fccols[(SubSet[j,i]/max(LogRatios,na.rm=T))*500+500], border=0)
                  }})
                text(0,0,"Log\nratios",cex=0.7)
                }
              }
            },height=600)
            
            incProgress(0.8, detail = paste("Plotting more results"))
            output$plotregdistr <- renderPlot({
              input$button
              qlim <- input$qval
              WhereRegs <-Qvalue
              WhereRegs[!is.null(WhereRegs)] <- 0
              WhereRegs[Qvalue<qlim] <- 1
              WhereRegs <- WhereRegs[,rep(0:4, NumCond-1)*(NumCond-1)+rep(1:(NumCond-1),each=5)]
              tcols <- rep(rainbow(NumCond-1),each=1)
             names(tcols) = rep(paste("A",1:(NumCond-1)),1)
              upset(as.data.frame(WhereRegs),nsets=ncol(WhereRegs),mainbar.y.label = "Significant features",order.by="degree",
                    decreasing=T,nintersects = NA,keep.order=T,sets=colnames(WhereRegs),text.scale=1.5, mb.ratio = c(0.55, 0.45),
                                   set.metadata = list(data = data.frame(set=colnames(WhereRegs),cols=paste("A",rep(1:(NumCond-1),each=5)),crab=1:ncol(WhereRegs)) , 
                                         plots = list(list(type = "matrix_rows",column = "cols", colors=tcols,alpha=0.5))))
            },height=600)
            
            output$plotreg <- renderPlot({
              input$button
              qlim <- input$qval
              par(mfrow=c(1,NumCond-1))
              tmpX <- 10^seq(log10(min(Qvalue,na.rm=T)),0.1,0.01)
              for (i in 1:(NumCond-1)) {
                plot(tmpX,rowSums(sapply((Qvalue[,i]),"<",tmpX),na.rm=T), main=paste("Comparison",i),xlab="q-value threshold",
                     ylab="Number significant",type="l",col="#3333FF",ylim=c(1,nrow(Qvalue)),log="xy",lwd=2)
                if (i==1)
                  legend("topleft",legend = testNames,col=c("#3333FF","#33FF33","#ff3333","#FFFF33","#FF33FF"),lwd=2)
                lines(tmpX,rowSums(sapply((Qvalue[,(NumCond-1)+i]),"<",tmpX),na.rm=T),col="#33FF33",lwd=2)
                lines(tmpX,rowSums(sapply((Qvalue[,(NumCond-1)*2+i]),"<",tmpX),na.rm=T),col="#FF3333",lwd=2)
                lines(tmpX,rowSums(sapply((Qvalue[,(NumCond-1)*3+i]),"<",tmpX),na.rm=T),col="#FFFF33",lwd=2)
                lines(tmpX,rowSums(sapply((Qvalue[,(NumCond-1)*4+i]),"<",tmpX),na.rm=T),col="#FF33FF",lwd=2)
                abline(v=qlim,col=2)
              }
              par(mfrow=c(1,1))
              
            },height=400)
            incProgress(0.9, detail = paste("Finishing"))
            
            output$downloadData <- downloadHandler(
              filename = function() {
                paste("Results", Sys.Date(), ".csv", sep="");
              },
              content = function(file) {
                write.csv(FullReg, file)
              })
            # output$downloadFigure <- downloadHandler(
            #   filename = function() {
            #     paste("Results", Sys.Date(), ".pdf", sep="");
            #   },
            #   content = function(file) {
            #     pdf(file,height=(NumCond-1)*4)
            #     print(dev.cur())
            #     replayPlot(pl)
            #     dev.off()
            #   })                
            
          })
          
          
        }
      })
    }
    #       output$messages <- renderText("finished")
    
  },height=heightSize)
})
