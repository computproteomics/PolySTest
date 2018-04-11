# Use hommel test (see https://stats.stackexchange.com/questions/76643/combining-p-values-from-different-statistical-tests-applied-on-the-same-data) for details and then smallest p-value for results

# library(shinyBS)
library(DT)
library(circlize)
library(UpSetR)
library(matrixStats)
library(qvalue)
library(heatmaply)
# library(d3heatmap)
library(gplots)
source("HelperFuncs.R")


options(shiny.maxRequestSize=2000*1024^2)

shinyServer(function(input, output,clientData,session) {
  dat <- FullReg <- NULL
  NumTests <- 6
  TestCols <- c("#33AAAA","#33AA33","#AA3333","#AA33AA","#AAAA33","#3333AA")
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
    updateNumericInput(session,"NumCond",value=4)
    updateNumericInput(session,"NumReps",value=3)
    updateNumericInput(session,"refCond",value=1)
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
    ColQuant <- input$ColQuant
    isPaired <- input$is_paired
    # print(dat)
    
    if (!is.null(dat)) {
      if (dat == "EXAMPLE")  {
        dat <- read.csv("LiverAllProteins.csv",row.names=1)
      } else {
        FullReg <<- NULL
        delim <- input$delimiter
        if (delim == "tab")
          delim <- "\t"
        if (input$row.names){
          dat <- read.csv(input$in_file$datapath,header=input$is_header,sep=delim,dec=input$digits)
          output$input_stats <- renderText("Duplicated feature names in first column! You can avoid them by not using 'Row names'")
          validate(need(sum(duplicated(dat[,1]),na.rm=T)==0,""))
          rownames(dat) <- dat[,1]
          dat <- dat[,2:ncol(dat)]
          # updateSliderInput(session,"QuantCol",max=ncol(dat)-2)
          
          # delete row with empty name
          dat <- dat[!rownames(dat)=="",]
          if (ColQuant > 2) {
            addInfo <- dat[,1:(input$ColQuant-2),drop=F]
            for (c in 1:ncol(addInfo))
              addInfo[,c] <- as.character(addInfo[,c])
            # print(head(addInfo))
            dat <- dat[,-(1:(ColQuant-2))]
          }
        } else {
          dat <- read.csv(input$in_file$datapath,header=input$is_header,sep=delim,dec=input$digits)
          # updateSliderInput(session,"QuantCol",max=ncol(dat))
          
          if (input$ColQuant > 1) {
            addInfo <- dat[,1:(ColQuant-1),drop=F]
            for (c in 1:ncol(addInfo))
              addInfo[,c] <- as.character(addInfo[,c])
            
            dat <- dat[,-(1:(ColQuant-1))]
          }
        }
      }
      print(ncol(dat))
      tncol <- 20
      if (!is.null(addInfo)) {
        tncol <- ncol(dat) + ncol(addInfo)
      } else {
        tncol <- ncol(dat)
      }
      
      updateNumericInput(session,"ColQuant",max=tncol)
      updateNumericInput(session,"NumCond",max=ncol(dat))
      updateNumericInput(session,"NumReps",max=ncol(dat))
      updateNumericInput(session,"refCond",max=input$NumCond)
      
      output$input_stats <- renderText(paste(ifelse(mode(as.matrix(dat))!="numeric","<b>Wrong file format /setup</b></br>",""),
                                             ifelse(ncol(dat) != NumReps*NumCond,"<b>Column number doesn't fit with number of replicates and conditions!</b><br/>",""),
                                             "Number of features: ",nrow(dat),
                                             "<br/>Number of data columns in file:", ncol(dat),
                                             "<br/>Percentage of missing values:",
                                             round(sum(is.na(dat))/nrow(dat)/ncol(dat)*100,digits = 2),"<br/>",
                                             paste("<br/>Comparison ",1:(NumCond-1),": Condition C",
                                                   (1:NumCond)[-(input$refCond)]," versus C",input$refCond,collapse=""),"<br/>"))
      
      ## Preview of input table
      Against<-seq(refCond,NumCond*NumReps,NumCond)
      RR<-1:(NumReps*NumCond)
      RR<-RR[-Against]
      RR<-rbind(RR,rep(Against,each=NumCond-1))
      compNames <- paste("C",RR[1,1:(NumCond-1)]," vs C",RR[2,1:(NumCond-1)],sep="")
      
      # Arrange table header
      sketch = htmltools::withTags(table(
        class = 'display',
        thead(
          tr(
            th('',style="text-align: center;"),
            if(!is.null(addInfo))
              th(colspan = ncol(addInfo), 'Metadata',style="text-align: center;border-left:thin solid;"),
            th(colspan = min(ncol(dat),NumReps*NumCond), 'Quantitative data columns',style="text-align: center;border-left:thin solid;")
          ),
          tr(
            th('',style="text-align: center;"),
            if(!is.null(addInfo))
              th(colspan = ncol(addInfo), '',style="text-align: center;border-left:thin solid;"),
            lapply(paste("C",rep(1:NumCond,NumReps)," Rep ",rep(1:NumReps,each=NumCond),sep="")[1:min(ncol(dat),NumCond*NumReps)],th,style="text-align: center;border-left:thin solid;")
          ),
          tr(
            th('Feature',style="text-align: center;"),
            if(!is.null(addInfo))
              lapply(colnames(addInfo),th,style="text-align: center;border-left:thin solid;"),
            lapply(colnames(dat)[1:min(ncol(dat),NumReps*NumCond)], th,style="text-align: center;border-left:thin solid;")
          )
        )
      ))
      if (!is.null(addInfo)) outTable <- cbind(addInfo,dat[,1:min(ncol(dat),NumReps*NumCond)])
      else outTable <- dat[,1:min(ncol(dat),NumReps*NumCond)]
      output$stat_table <- DT::renderDataTable(DT::datatable(outTable,
                                                             filter = list(position = 'top', clear = FALSE),colnames = c('model' = 1),
                                                             options = list(scrollX = TRUE,dom = 'Blfrtip',
                                                                            columnDefs = list(list(width = '20%', targets = colnames(data))),
                                                                            autoWidth=T),
                                                             class="compact",
                                                             container=sketch))
      
      #print(head(dat))
      
      
      isolate({
        if (input$button == 0)
          return()
        if (ncol(dat) == NumReps*NumCond) {
          withProgress(message="Calculating ...", min=0,max=1, {
            
            #           output$messages <- renderText("running")
            incProgress(0.1, detail = paste("Running statistical tests"))
            dat [!is.finite(as.matrix(dat))] <- NA
            dat <<- dat
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
            
            
            setProgress(0.5, detail = paste("Preparing data"))
            
            LogRatios <- qvalues$lratios
            Pvalue <- cbind(qvalues$plvalues, qvalues$pRPvalues, qvalues$pPermutvalues, MissingStats$pNAvalues, qvalues$ptvalues)
            Qvalue <- cbind(qvalues$qlvalues, qvalues$qRPvalues, qvalues$qPermutvalues, MissingStats$qNAvalues,qvalues$qtvalues)
            Qvalue <- cbind(UnifyQvals(Qvalue,NumCond,NumTests),Qvalue)
            # WhereReg <- cbind(qvalues$qtvalues<qlim, qvalues$qlvalues<qlim, qvalues$qRPvalues<qlim, qvalues$qPermutvalues<qlim, MissingStats$qNAvalues<qlim)
            
            testNames <- c("limma","rank products","Permutation test","NA test","t-test")
            colnames(LogRatios) <- paste("log-ratios",compNames)
            colnames(Pvalue) <- paste("p-values",rep(testNames,each=NumCond-1),rep(compNames,length(testNames)))
            testNames2 <- c("unified tests",testNames)
            names(TestCols) <- testNames2
            updateCheckboxGroupInput(session,"selTests",choices = testNames2, selected = "unified tests")
            updateCheckboxGroupInput(session,"selComps",choices = compNames, selected=compNames[1])
            colnames(Qvalue) <- paste("q-values",rep(testNames2,each=NumCond-1),rep(compNames,length(testNames2)))
            # colnames(WhereReg) <- paste("Differentially regulated",rep(testNames,each=NumCond-1),rep(compNames,length(testNames)))
            
            # print(cor(log10(Pvalue),use="na.or.complete"))
            if (!is.null(addInfo))
              FullReg <- cbind(addInfo[rownames(LogRatios),],LogRatios, Qvalue)#, WhereReg)
            else 
              FullReg <- cbind(LogRatios, Qvalue)#, WhereReg)
            
            # Calculate best fc and qlim combination for all tests but the t-test
            setProgress(0.7, detail = paste("Calculating favorable q-value and fc thresholds"))
            tcomb <- FindFCandQlim(Qvalue, LogRatios)
            setProgress(0.9, detail = paste("Creating figures"))
            
            print(tcomb)
            
            # Set fold-change slider range
            updateSliderInput(session,"fcval",min=round(min(LogRatios,na.rm=T),1),
                              max=round(max(LogRatios,na.rm=T),1), value=c(-tcomb[1],tcomb[1]))
            
            # Set q-value threshold
            updateNumericInput(session,"qval",value=tcomb[2])
            
            output$table_stats <- renderText(paste("Number selected features:",length(input$stat_table_rows_selected)))
            
            
            # Arrange table header
            sketch = htmltools::withTags(table(
              class = 'display',
              thead(
                tr(
                  th('',style="text-align: center;"),
                  if(!is.null(addInfo))
                    th(colspan = ncol(addInfo), 'Metadata',style="text-align: center;border-left:thin solid;"),
                  th(colspan = NumCond-1, 'Log-ratios',style="text-align: center;border-left:thin solid;"),
                  th(colspan = (NumCond-1)*NumTests, 'q-values',style="text-align: center;border-left:thin solid;")
                ),
                tr(
                  th('',style="text-align: center;"),
                  if(!is.null(addInfo))
                    th(colspan = ncol(addInfo), ''),
                  th(colspan = NumCond-1, '',style="text-align: center;border-left:thin solid;border-left:thin solid;"),
                  th(colspan = (NumCond-1), 'unified',style="text-align: center;border-left:thin solid;color: #AA3333;"),
                  th(colspan = (NumCond-1), 'limma',style="text-align: center;border-left:thin solid;"),
                  th(colspan = (NumCond-1), 'rank products',style="text-align: center;border-left:thin solid;"),
                  th(colspan = (NumCond-1), 'permutation test',style="text-align: center;border-left:thin solid;"),
                  th(colspan = (NumCond-1), 'NA test',style="text-align: center;border-left:thin solid;"),
                  th(colspan = (NumCond-1), 't-test',style="text-align: center;border-left:thin solid;")
                ),
                tr(
                  th('Feature',style="text-align: center;"),
                  if(!is.null(addInfo))
                    lapply(colnames(addInfo),th,style="text-align: center;border-left:thin solid;"),
                  lapply(rep(compNames,NumTests+1), th,style="text-align: center;border-left:thin solid;")
                )
              )
            ))
            output$stat_table <- DT::renderDataTable(DT::datatable(FullReg,
                                                                   filter = list(position = 'top', clear = FALSE),colnames = c('model' = 1),
                                                                   options = list(scrollX = TRUE,dom = 'Blfrtip',
                                                                                  columnDefs = list(list(width = '20%', targets = colnames(FullReg))),
                                                                                  autoWidth=T,lengthMenu = c(5, 10, 50,100),
                                                                                  buttons = list('colvis', 'copy', 'print')),
                                                                   extensions=c("Buttons","FixedColumns"),class="compact",
                                                                   container=sketch) %>% 
                                                       formatSignif(grep("log-ratios",colnames(FullReg),value=T),digits=2) %>%
                                                       formatSignif(grep("q-values",colnames(FullReg),value=T),digits=2))
            
            
            SubSetLR <- SubSetQval <- FCRegs <- NULL
            qlim <- 0.01
            fclim <- c(-1,1)
            
            proxy = dataTableProxy('stat_table')
            
            observeEvent(input$resetSelection, {
              proxy %>% DT::selectRows(NULL)
            })
            
            observeEvent(input$allLimsSelection, {
              selCols <- paste("q-values",apply(expand.grid(input$selTests, input$selComps), 1, paste, collapse=" "))
              print(colnames(FCRegs))
              selFeat <- which(rowSums(FCRegs[,selCols,drop=F]<input$qval)>0)
              proxy %>% DT::selectRows(NULL)
              if (length(selFeat) > 0)
                proxy %>% DT::selectRows(as.numeric(selFeat))
            })
            
            observeEvent(input$allPageSelection, {
              proxy %>% DT::selectRows(NULL)
              proxy %>% DT::selectRows(input$stat_table_rows_current)
            })
            
            observeEvent(input$allSelection, {
              proxy %>% DT::selectRows(NULL)
              proxy %>% DT::selectRows(input$stat_table_rows_all)
            })
            
            observeEvent(input$plotregdistr_click, {
              print(input$plotregdistr_click$x)
              print(input$plotregdistr_click$y)
              
            })
            
            ## Delay reaction to selecting rows in data table
            triggerUpdate <- debounce(reactive(input$stat_table_rows_selected),1000)
        
            observe({
              input$button
              input$stat_table
              sel_prots <- triggerUpdate()
              qlim <- input$qval
              fclim <- input$fcval
              # print(input$stat_table_rows_selected)
              # Selecting only features from selected tests and conditions
              # print(input$selComps)
              # print(input$selTests)
              SubSetLR <<- LogRatios[sel_prots,,drop=F]
              SubSetQval <<- Qvalue[sel_prots,,drop=F]
              SubSetLR <<- SubSetLR[order(rowMins(SubSetQval[,1:(NumCond-1),drop=F],na.rm=T)),,drop=F]
              SubSetQval <<- SubSetQval[order(rowMins(SubSetQval[,1:(NumCond-1),drop=F],na.rm=T)),,drop=F]
              ## Same as Qvalue but with NAs and corresponding fold-changes filtered and set to 1
              FCRegs <<- Qvalue
              for (t in 1:NumTests) {
                tsign <- Qvalue[,(t-1)*(NumCond-1)+(1:(NumCond-1)),drop=F] 
                tsign[is.na(tsign)] <- 1
                tsign[LogRatios>fclim[1] & LogRatios<fclim[2]] <- 1
                FCRegs[,(t-1)*(NumCond-1)+(1:(NumCond-1))] <<- tsign
              }
            })
            
            
            
            incProgress(0.7, detail = paste("Plotting first results"))
            
            
            ## Plotting DRFs vs q-value thresholds + volcano plots + p-val vs q-vals
            plotPvalueDistr <- function() {
              par(mfrow=c(NumCond-1,length(testNames)))
              for (i in 1:(NumCond-1)) {
                for (j in 2:NumTests) {
                  hist(Pvalue[,(NumCond-1)*(j-2)+i],100, main=testNames2[j],sub=compNames[i],col=TestCols[j],xlab="p-value",border=NA)
                }
              }
            }
            plotPvalueDistr()
            output$downloadPvalueDistrPdf <- downloadHandler(
              filename = function() {
                paste("PvalueDistrPlots", Sys.Date(), ".pdf", sep="");
              },
              content = function(file) {
                pdf(file,height=12,width=20)
                plotPvalueDistr()
                dev.off()  
              })
            
            par(mfrow=c(1,1))
            # volcano plots
            output$plotvolc <- renderPlot({
              input$button
              input$stat_table
              sel_prots <- triggerUpdate()
              qlim <- input$qval
              fclim <- input$fcval
              plotVolcano <- function() {
                par(mfrow=c(NumCond-1,NumTests))
                for (i in 1:(NumCond-1)) {
                  for (j in 1:NumTests) {
                    plot(LogRatios[,i],-log10(Qvalue[,(NumCond-1)*(j-1)+i]), main=testNames2[j],sub=compNames[i],xlab="log fold-change",ylab="-log10(q)",
                         cex=colSelected(0.5,nrow(Qvalue),sel_prots,1),
                         col=colSelected(adjustcolor(TestCols[j],alpha.f=0.3),nrow(Qvalue),sel_prots,"#FF9933"),pch=16,ylim=-log10(c(1,min(Qvalue,na.rm=T))))
                    abline(h=-log10(qlim),col="#AA3333",lwd=2)
                    abline(v=fclim,col="#AA3333",lwd=2)
                  }
                }
              }
              plotVolcano()
              par(mfrow=c(1,1))
              output$downloadVolcanoPdf <- downloadHandler(
                filename = function() {
                  paste("VolcanoPlots", Sys.Date(), ".pdf", sep="");
                },
                content = function(file) {
                  pdf(file,height=12,width=20)
                  plotVolcano()
                  dev.off()  
                })
              
              
            },height=heightSize)

            output$plotexpression <- renderPlot({
              input$button
              input$stat_table
                                          
              triggerUpdate()
              qlim <- input$qval
              fclim <- input$fcval
              #input$stat_table_rows_selected
              # print(head(SubSetLR))
              if (length(SubSetLR)> 0) {
                # CI plots of max 30 features
                SubSet <- SubSetLR[1:min(nrow(SubSetLR),30),,drop=F]
                indices <- rownames(SubSet)
                tdat <- as.matrix(dat[rownames(SubSet), (rep(1:NumReps,NumCond)-1)*NumCond+rep(1:NumCond,each=NumReps),drop=F])
                MeanSet <- SDSet <-  matrix(NA,nrow=nrow(tdat),ncol=NumCond,dimnames = 
                                              list(x=rownames(tdat),y=paste("Condition",1:NumCond)))
                for (c in 1:NumCond) {
                  MeanSet[,c] <- rowMeans(tdat[,1:NumReps + (c-1)*NumReps,drop=F],na.rm=T)
                  SDSet[,c] <- rowSds(tdat[,1:NumReps + (c-1)*NumReps,drop=F],na.rm=T)
                }
                # par(mfrow=c(1,3))
                plotExpression <- function() {
                  layout(t(c(1,1,2,2,3,3)))
                  plot(0,0,type="n",bty="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA)
                  # print(colnames(SubSet))
                  legend("topright",col=rainbow(nrow(SubSet),alpha = 0.8,s=0.7),legend=rownames(SubSet),lwd=3)
                  # plotCI(1:(NumCond-1)+runif(1,-0.1,0.1),LogRatios[as.vector(indices)[1],1:(NumCond-1),drop=F],pch=16,
                  #        xlab="Conditions",xlim=c(0.7,(NumCond)-0.7),
                  #        ylab="log-ratios",col=rainbow(nrow(SubSet),alpha = 0.8,s=0.7)[1],
                  #        uiw=qvalues$Sds[input$stat_table_rows_selected][1],type="b",barcol="#000000FF",
                  #        ylim=range(LogRatios,na.rm=T),xaxt="none",lwd=1.5)
                  plotCI(1:(NumCond)+runif(1,-0.1,0.1),MeanSet[1,],pch=16,
                         xlab="Conditions",xlim=c(0.7,(NumCond+1)-0.7),
                         ylab="expression values",col=rainbow(nrow(MeanSet),alpha = 0.8,s=0.7)[1],
                         uiw=SDSet[1,],type="b",barcol="#000000FF",
                         ylim=range(tdat,na.rm=T),xaxt="none",lwd=1.5)
                  axis(1,at=1:(NumCond),labels = colnames(MeanSet))
                  # abline(h=fclim)
                  if (nrow(MeanSet)>1) {
                    for (i in 2:nrow(MeanSet)){
                      # plotCI(1:(NumCond-1)+runif(1,-0.1,0.1),LogRatios[as.vector(indices[i]),1:(NumCond-1),drop=F],
                      #        add = T,pch=16,col=rainbow(nrow(SubSet))[i],
                      #        uiw=qvalues$Sds[input$stat_table_rows_selected][i],type="b",barcol="#000000FF",lwd=1.5) 
                      plotCI(1:(NumCond)+runif(1,-0.1,0.1),MeanSet[i,],
                             add = T,pch=16,col=rainbow(nrow(MeanSet))[i],
                             uiw=SDSet[i,],type="b",barcol="#000000FF",lwd=1.5) 
                      
                    }
                  }
                  
                  # Missing values
                  # find descent visualization 
                  
                  # Compare q-values for each protein and method
                  # validate(
                  #   need(length(SubSet)>0, "Please select features from the data table")
                  # )
                  
                  if (length(SubSet)> 0) {
                    
                    par(mar=rep(0,4))
                    circos.clear()
                    circos.par(cell.padding=c(0,0,0,0),canvas.xlim=c(-1.5,1.5),canvas.ylim=c(-1.5,1.5),
                               track.margin=c(0,0.02),start.degree=90,gap.degree=4)
                    circos.initialize(1:(NumCond-1),xlim=c(0,1))
                    for (t in 1:NumTests) {
                      # print(tsign)
                      nfeat <- min(nrow(SubSet),30)
                      cols <- rainbow(nfeat,alpha = 0.8,s=0.7)
                      tsign <- FCRegs[indices,(t-1)*(NumCond-1)+(1:(NumCond-1)),drop=F]<qlim
                      circos.trackPlotRegion(ylim=c(-3,2),track.height=1/12, bg.border="#777777",
                                             panel.fun = function(x, y) {
                                               name = get.cell.meta.data("sector.index")
                                               i = get.cell.meta.data("sector.numeric.index")
                                               xlim = get.cell.meta.data("xlim")
                                               ylim = get.cell.meta.data("ylim")
                                               xdiff <- (xlim[2]-xlim[1])/nfeat
                                               if(t == 1) {
                                                 circos.text(mean(xlim), max(ylim)+30, compNames[i], facing = "inside", niceFacing = TRUE,cex = 1,font=2)
                                                 circos.axis("top", labels = rownames(SubSetLR),major.at=seq(1/(nfeat*2),1-1/(nfeat*2),length=nfeat),minor.ticks=0,
                                                             labels.cex = 0.8,labels.facing = "reverse.clockwise")
                                               }
                                               for (j in which(tsign[,i])) {
                                                 circos.rect(xleft=xlim[1]+(j-1)*xdiff, ybottom=ylim[1],
                                                             xright=xlim[2]-(nfeat-j)*xdiff, ytop=ylim[2],
                                                             col = cols[j], border=cols[j])
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
                    mtext(paste("Track ",1:NumTests,": ",testNames2,sep="",collapse="\n"),
                          side=1,outer=T,adj=1,line=-1,cex=0.7)
                  }
                }
                plotExpression()
                
                output$downloadExprPdf <- downloadHandler(
                  filename = function() {
                    paste("SelExpressionProfiles", Sys.Date(), ".pdf", sep="");
                  },
                  content = function(file) {
                    pdf(file,height=8,width=12)
                    plotExpression()
                    dev.off()  
                  })
                
              }
            },height=400)
            
            output$plotheatmap <- renderPlotly({
              # d3heatmap(SubSetLR)
              input$button
              input$stat_table
              triggerUpdate()
              qlim <- input$qval
              fclim <- input$fcval
              print((SubSetLR))
              # SubSetLR <- SubSetLR[rowSums(!is.na(SubSetLR))>1,,drop=F]
              p <- plotly_empty()
              if (!is.null(SubSetLR)) {
                if (length(SubSetLR)> 0 & nrow(SubSetLR)>1) {
                  print("running heatmap")
                  withProgress(message="Creating heatmap ...", min=0,max=1, {
                    setProgress(0.5)
                    # print(SubSetLR)
                    tdat <- dat[rownames(SubSetLR), (rep(1:NumReps,NumCond)-1)*NumCond+rep(1:NumCond,each=NumReps)]
                    # remove data rows with more than 45% missing values
                    to_remove <- which(rowSums(is.na(tdat)) > ncol(tdat)*0.45)
                    tqvals <- Qvalue[rownames(SubSetLR),1:(NumCond-1)]
                    if (length(to_remove)>0) {
                      tqvals <- tqvals[-to_remove,]
                      tdat <- tdat[-to_remove,]
                    }
                    # setting colors of p-values
                    pcols <- rev(c(0.001, 0.01, 0.05, 1))
                    ttt <- tqvals
                    for (c in pcols) {
                      ttt[tqvals <= c] <- c 
                    }
                    tqvals <- data.frame(ttt)
                    for (c in 1:ncol(tqvals))
                     tqvals[,c] <- paste("<",as.character(ttt[,c],pcols),sep="")

                    p <- heatmaply(tdat[order(rownames(tdat)),],Colv=F,scale = "none",trace="none",cexRow=0.7,plot_method="plotly", 
                                   RowSideColors = tqvals, row_side_palette = grey.colors)
                    # p <- heatmaply(SubSetLR,scale = "none",trace="none",cexRow=0.7)
                    # heatmap.2(SubSetLR,col=bluered,cexCol = 0.7,srtCol=45,scale="none",trace="none",cexRow=0.7)
                  })
                }
              }
              p
              
            })#,height=800)
            
            incProgress(0.8, detail = paste("Plotting more results"))
            output$plotregdistr <- renderPlot({
              qlim <- input$qval
              input$fcval
              input$button
              # print(head(FCRegs))
              WhereRegs <- FCRegs[,rep(0:(NumTests-2), NumCond-1)*(NumCond-1)+rep(1:(NumCond-1),each=NumTests-1),drop=F]<qlim
              # print(head(WhereRegs))
              WhereRegs[WhereRegs] <- 1
              deleted_cols <- which(colSums(WhereRegs,na.rm=T)==0)
              # print(deleted_cols)
              
              tcolnames <- paste("A",rep(1:(NumCond-1),each=NumTests-1))
              if (length(deleted_cols) > 0) {
                tcolnames <- tcolnames[-deleted_cols]
                WhereRegs <- WhereRegs[,-deleted_cols,drop=F]
              }
              tcols <- rep(rainbow(NumCond-1),each=1)
              names(tcols) = rep(paste("A",1:(NumCond-1)),1)
              # print(head(WhereRegs))
              plotUpset <- function () {
                upset(as.data.frame(WhereRegs),nsets=ncol(WhereRegs),mainbar.y.label = "Significant features",order.by="degree",
                      decreasing=T,nintersects = NA,keep.order=T,sets=colnames(WhereRegs),text.scale=1.5, mb.ratio = c(0.55, 0.45),
                      set.metadata = list(data = data.frame(set=colnames(WhereRegs),cols=tcolnames,crab=1:ncol(WhereRegs)), 
                                          plots = list(list(type = "matrix_rows",column = "cols", colors=tcols,alpha=0.5))))
              }
              plotUpset()
              output$downloadUpSetPdf <- downloadHandler(
                filename = function() {
                  paste("UpSetProfiles", Sys.Date(), ".pdf", sep="");
                },
                content = function(file) {
                  pdf(file,height=8,width=8)
                  plotUpset()
                  dev.off()  
                })
              
            },height=600)
            
            output$plotreg <- renderPlot({
              qlim <- input$qval
              input$fcval
              input$button
              par(mfrow=c(1,NumCond-1))
              tmpX <- 10^seq(log10(min(Qvalue,na.rm=T)),0.1,0.01)
              plotRegDistr <- function() {
                for (i in 1:(NumCond-1)) {
                  plot(tmpX,rowSums(sapply((FCRegs[,i]),"<",tmpX),na.rm=T), main=paste("Comparison",i),xlab="q-value threshold",
                       ylab="Number significant",type="l",col=TestCols[1],ylim=c(1,nrow(Qvalue)),log="xy",lwd=2)
                  if (i==1)
                    legend("topleft",legend = testNames2,col=TestCols,lwd=2)
                  lines(tmpX,rowSums(sapply((FCRegs[,(NumCond-1)*1+i]),"<",tmpX),na.rm=T),col=TestCols[2],lwd=2)
                  lines(tmpX,rowSums(sapply((FCRegs[,(NumCond-1)*2+i]),"<",tmpX),na.rm=T),col=TestCols[3],lwd=2)
                  lines(tmpX,rowSums(sapply((FCRegs[,(NumCond-1)*3+i]),"<",tmpX),na.rm=T),col=TestCols[4],lwd=2)
                  lines(tmpX,rowSums(sapply((FCRegs[,(NumCond-1)*4+i]),"<",tmpX),na.rm=T),col=TestCols[5],lwd=2)
                  lines(tmpX,rowSums(sapply((FCRegs[,(NumCond-1)*5+i]),"<",tmpX),na.rm=T),col=TestCols[6],lwd=2)
                  abline(v=qlim,col=2)
                }
              }
              plotRegDistr()
              par(mfrow=c(1,1))
              output$downloadRegDistrPdf <- downloadHandler(
                filename = function() {
                  paste("RegDistrPlots", Sys.Date(), ".pdf", sep="");
                },
                content = function(file) {
                  pdf(file,height=8,width=8)
                  plotRegDistr()
                  dev.off()  
                })
              
              
            },height=400)
            setProgress(0.9, detail = paste("Finishing"))
            
            output$downloadData <- downloadHandler(
              filename = function() {
                paste("Results", Sys.Date(), ".csv", sep="");
              },
              content = function(file) {
                write.csv(cbind(FullReg,Selected=(1:nrow(FullReg) %in% input$stat_table_rows_selected),dat), file)
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