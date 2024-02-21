# Use hommel test (see https://stats.stackexchange.com/questions/76643/combining-p-values-from-different-statistical-tests-applied-on-the-same-data)
# for details and then smallest p-value for results

library(DT)
library(circlize)
library(UpSetR)
library(matrixStats)
library(qvalue)
library(heatmaply)
library(SummarizedExperiment)
library(gplots)
library(jsonlite)
library(PolySTest)

validate <- shiny::validate

options(shiny.maxRequestSize=2000*1024^2)
options(java.parameters="-Xss2560k")
# for debugging:
# options(shiny.reactlog=TRUE)

shinyServer(function(input, output, clientData, session) {

  # Reactive values and global variables
  currdata <- reactiveVal(NULL)  # Stores the current data content from file
  currse <- reactiveVal(NULL) # Stores data in SummarizedExperiment
  FullReg <- reactiveVal(NULL)  # Output table
  obsadd <- NULL  # Placeholder for adding new comparison UI elements
  actFileName <- NULL  # Placeholder for the actual file name for input
  extdatatable <- reactiveVal(NULL)  # Stores external data table submitted to PolySTest via JSON message
  Comps <- reactiveValues(ind=vector(), num=NULL, RR=NULL, compNames=NULL, allComps=NULL)  # Stores comparison information
  currButton <- 0  # Current button state of input$button
  addCompButton <- 0  # Check whether new comparison button needs to be added
  NumTests <- 6  # Number of different statistical tests including the combined PolySTest
  TestCols <- c("#33AAAA","#33AA33","#AA3333","#AA33AA","#AAAA33","#3333AA")  # Colors for tests
  testNames <- c("limma","Miss test","rank products","permutation test","t-test")
  testNames2 <- c("PolySTest",testNames)
  names(TestCols) <- testNames2
  addInfo <- NULL  # Additional information from the original table to be re-added after analysis
  proxy <- NULL # Access to data table
  heightSize <- reactiveVal(200)  # Height of the data table



  # Render text outputs
  output$messages <- renderText("")
  output$thead <- renderText("Multiple statistical tests for the detection
                             for differentially regulated features")
  output$description <- renderText("The methods are designed to carry out
                                   statistical tests on data with few
                                   replicates and high amounts of missing values.")
  output$description2 <- renderText("We recommend this method for datasets with a minimum of 1000 features (rows of the table) and at least 3 replicates. The method is based on ratios between features within the replicate, i.e. the tests are carried out on paired tests.")
  output$description3 <- renderText("For details, see <i>manuscript in preparation</i>")
  output$description4 <- renderText("Data input: Data table of <b>log-transformed</b> (quantitative) omics data (csv file) with separate columns for feature information and quantification values. If row names is set, then the first column of the csv file provides unique features names (no duplicates). Columns with quantifications follow after all other columns and start with columns number <i>First column for quantification</i>. The quantification values should be log-transformed intensity/abundance values (not ratios) which have been normalized to be comparable. Be aware that not normalized or wrongly normalized data might increase the number of false positives. The order of the columns is required to either be A1, A2, A3, ..., B1, B2, B3, ... or A1, B1, C1, ..., where 1,2,3 ... are the conditions and A,B,... denote replicates. We require identical number of replicates per condtion. In the case of slight differing replicate numbers, add empty columns. The tests check for differentially regulated features versus the \"reference\" condition. For each comparison (log-ratio), we provide various visualizations for further investigation.")

  # Download data handler
  output$downloadData <- downloadHandler(
    filename = function() {
      validate(
        need(F, "No data")
      )
      paste("Results", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(FullReg(), file)
    }
  )

  # Observe event for selecting example data
  observeEvent(input$example, {
    currdata("EXAMPLE")
  })


  # load external data if submitted
  observe({
    if (!is.null(input$extdata)) {
      isolate({
        jsonmessage <- fromJSON(input$extdata)
        print(js$send_results)
        # Loading parameters
        NumCond <- jsonmessage[["numcond"]]
        NumReps <- jsonmessage[["numrep"]]
        isPaired <- jsonmessage[["paired"]]
        isGrouped <- jsonmessage[["grouped"]]
        ColQuant <- jsonmessage[["firstquantcol"]]
        # reading data matrix
        expr_matr <- jsonmessage[["expr_matrix"]]
        print(names(expr_matr))
        output$fileInText <- renderText({
          validate(need(!is.null(expr_matr), "Uploaded data empty"))
          validate(need(length(expr_matr) > 1, "Uploaded data does not contain multiple columns"))
          validate(need(sum(duplicated(expr_matr[[1]]), na.rm = T) == 0, "Duplicated feature names in first column!"))
          tdat <- expr_matr[[1]]
          for (i in 2:length(expr_matr)) {
            validate(need(length(expr_matr[[i]]) == length(expr_matr[[1]]),
                          paste("Wrong array length of sample", names(expr_matr)[i])))
            if (i < ColQuant) {
              tdat <- data.frame(tdat, expr_matr[[i]])
            } else {
              tdat <- data.frame(tdat, as.numeric(expr_matr[[i]]))
            }
          }
          colnames(tdat) <- names(expr_matr)
          updateNumericInput(session, "NumCond", value = NumCond)
          updateNumericInput(session, "NumReps", value = NumReps)
          updateCheckboxInput(session, "is_paired", value = isPaired)
          updateCheckboxInput(session, "qcol_order", value = isGrouped)
          updateCheckboxInput(session, "ColQuant", value = ColQuant)
          extdatatable(tdat)
          print("Loaded external data")
          return(paste("Loaded external data"))
        })

      })
    }
  })


  # adding default UI elements for comparisons
  addCompUIs <- function(el, conditions) {
    print("addCompUIs")
    div(style="padding-right: 10px; padding-left: 0px;",id=paste("selall_",el,sep=""),
        fluidRow(
          column(12,align="left",style="padding:0px;",div(p(paste("Comparison ", el, ":",sep="")),h6(style = "display:inline;", icon("question-circle"))),
                 id=paste("selt_",el,sep=""),style="padding:0px;")
        ),
        fluidRow(
          column(4,align="center",style="padding:0px;",selectInput(paste("sels_",el,sep=""),label=NULL,
                                                                   choices = conditions,selected = conditions[el+1])),
          column(2,h5("vs")),
          column(4,align="center",style="padding:0px;",selectInput(paste("selr_",el,sep=""),
                                                                   label=NULL,
                                                                   choices = conditions,selected = conditions[1])),
          column(2,align="center",style="padding:0px;",actionButton(paste("selb_",el,sep=""),label=NULL,icon =icon("trash")))
        ))
  }

  observe({
    # input$example
    print("observe, read")
    output$messages <- renderText("")
    dat <- currdata()
    if (!is.null(input$in_file)) {
      dat <- input$in_file
      print("reading file")
      delim <- input$delimiter
      if (delim == "tab")
        delim <- "\t"
      if (length(dat) == 4) {
        actFileName <<- input$in_file$datapath
      }
      dat <- read.csv(actFileName,header=input$is_header,sep=delim,dec=input$digits,stringsAsFactors = F)
    } else if (!is.null(extdatatable())) {
      dat <- extdatatable()
      print(head(data))
    } else if (!is.null(actFileName)) {
      print("reading file")
      delim <- input$delimiter
      if (delim == "tab")
        delim <- "\t"
      dat <- read.csv(actFileName,header=input$is_header,sep=delim,dec=input$digits,stringsAsFactors = F)
    }
    NumReps <- input$NumReps
    NumCond <- input$NumCond
    if(is.na(NumReps)) NumReps <- 1
    if(is.na(NumCond)) NumCond <- 1

    ColQuant <- input$ColQuant
    isPaired <- input$is_paired
    print(head(dat,n=1))
    validate(need(!is.na(input$ColQuant),"Change 'First column for quantification'"))
    if (!is.null(dat)) {
      ndatcol <- 0
      if (unlist(dat)[1] == "EXAMPLE")  {
        actFileName <<- system.file("extdata", "LiverAllProteins.csv", package = "PolySTest")
        dat <- read.csv(actFileName,row.names=1)
        # limiting to only 500 rows
        dat <- dat[1:500, ]
        updateNumericInput(session,"NumCond",value=4)
        NumCond <- 4
        updateNumericInput(session,"NumReps",value=3)
        NumReps <- 3
        updateCheckboxInput(session,"qcol_order",value=T)
        updateCollapse(session,"Input",open = "Statistical testing",close="Data input")
        ndatcol <- 12
      } else  {
        FullReg(NULL)
        if (input$row.names){
          output$input_stats <- renderText("Duplicated feature names in first column! You can avoid them by not using the option 'Row names'")
          validate(need(sum(duplicated(dat[,1]),na.rm=T)==0,""))
          rownames(dat) <- dat[,1]
          dat <- dat[,2:ncol(dat)]
          # updateSliderInput(session,"QuantCol",max=ncol(dat)-2)
          # delete row with empty name
          dat <- dat[!rownames(dat)=="",]
          if (input$ColQuant > 2) {
            validate(need(input$ColQuant-2 <= ncol(dat), "To large 'First column for quantification'!"))
            addInfo <<- dat[,1:(input$ColQuant-2),drop=F]
            for (c in 1:ncol(addInfo))
              addInfo[,c] <- as.character(addInfo[,c])
            # print(head(addInfo))
            dat <- dat[,-(1:(input$ColQuant-2))]
          } else {
            addInfo <- NULL
          }

        } else {
          # updateSliderInput(session,"QuantCol",max=ncol(dat))
          rownames(dat) <- paste("feature",1:nrow(dat))

          if (input$ColQuant > 1) {
            validate(need(input$ColQuant-1 <= ncol(dat), "To large 'First column for quantification'!"))
            addInfo <<- dat[,1:(input$ColQuant-1),drop=F]
            for (c in 1:ncol(addInfo))
              addInfo[,c] <- as.character(addInfo[,c])

            dat <- dat[,-(1:(input$ColQuant-1))]
          }
        }

        ndatcol <- ncol(dat)
        if (!input$qcol_order) {
          print("reorder columns")
          validate(need(ncol(dat)>=NumCond*NumReps, "Not enough quantitative columns"))
          act_cols <- rep(0:(NumCond-1),NumReps)*NumReps+rep(1:(NumReps), each=NumCond)
          dat <- dat[,act_cols]
        }

        tncol <- 20
        if (!is.null(addInfo)) {
          tncol <- ncol(dat) + ncol(addInfo)
        } else {
          tncol <- ncol(dat)
        }
        updateNumericInput(session,"ColQuant",max=tncol)
        updateNumericInput(session,"NumCond",max=ncol(dat))
        updateNumericInput(session,"NumReps",max=ncol(dat))
        updateCollapse(session,"Input",open = "Data layout",close="Data input")

      }
      conditions <- paste("C",1:NumCond,sep="")

      for (i in 0:100) {
        removeUI(selector=paste("#selall_",i,sep=""),immediate=T)
      }
      Comps$num <- length(conditions)-1
      Comps$ind <- 1:(length(conditions)-1)
      print("Conditions:")
      for (el in (length(conditions)-1):1) {
        insertUI("#stat_comparisons","afterEnd",ui=tagList(addCompUIs(el,conditions)), immediate=T)
        addTooltip(session,id=paste("selt_",el,sep=""),title="Condition and reference condition which will be compared (taking log-ratios)",trigger="hover")
        addTooltip(session, paste("selall_",el,sep=""), "")
      }

      # add new UI elements for comparison
      if (!is.null(obsadd)) {
        print("obsadd")
        obsadd$destroy()
      }


      obsadd <<- observeEvent(input$addComp, {
        print("addComp")
        if (input$addComp == addCompButton)
          return()
        addCompButton <<- input$addComp


        ind <- 1
        if (length(Comps$ind)>0)
          ind <- min((1:100)[-Comps$ind])
        insertUI("#addComp","beforeBegin",ui=tagList(addCompUIs(ind,conditions)), immediate=T)
        addTooltip(session,id=paste("selt_",ind,sep=""),title="Condition and
                   reference condition which will be compared
                   (taking log-ratios)",trigger="hover")
        Comps$num <- Comps$num + 1
        Comps$ind <- c(Comps$ind, ind)

      })

      # remove UI elements
      lapply(0:100, function(i) {
        obs <- observeEvent(input[[paste("selb_",i,sep="")]],{
          isolate({
            # output[[paste("div:has(> #","selr_",i,")",sep="")]] <- NULL
            # output[[paste("div:has(> #","sels_",i,")",sep="")]] <- NULL
            # removeUI(selector=paste("div:has(> #","selr_",i,")",sep=""))
            # removeUI(selector=paste("div:has(> #","selb_",i,")",sep=""))
            # removeUI(selector=paste("div:has(> #","selt_",i,")",sep=""))
            print(paste("Remove sellall_",i))
            removeUI(selector=paste("#","selall_",i,sep=""),immediate=T)
            if (sum(Comps$ind == i) > 0) {
              Comps$num <- Comps$num - 1
              Comps$ind <- Comps$ind[-which(Comps$ind == i)]
            }
          })
        })

      })


      if(ncol(dat) == NumReps*NumCond) {
        updateCollapse(session,"Input",open="Statistical tests")
      } else {
        updateCollapse(session,"Input",close="Statistical tests")
      }

      ## Preview of input table

      isolate({
        output$input_stats <- renderText({
          print("table calc")
          allComps <- NULL
          for(cmp in Comps$ind) {
            if (!is.null(input[[paste("selr_",cmp,sep="")]]) & !is.null(input[[paste("sels_",cmp,sep="")]]))
              allComps <- rbind(allComps, c(input[[paste("selr_",cmp,sep="")]], input[[paste("sels_",cmp,sep="")]]))
          }
          allComps <- unique(allComps)
          allComps <- allComps[allComps[,1,drop=F] != allComps[,2,drop=F], , drop=F]
          compNames <- NULL
          ncomps <- nrow(allComps)

          if(!is.null(allComps)) {
            valComps <- matrix(NA,nrow=ncomps, ncol=ncol(allComps))
            for (i in 1:length(allComps)) valComps[i] <- as.numeric(sub("C","",allComps[i]))
            RR <- matrix(NA,ncol = ncomps*NumReps, nrow=2)
            for (j in 1:ncomps) {
              compCond <- valComps[j,2]
              refCond <- valComps[j,1]
              RR[1,seq(j,ncomps*NumReps,ncomps)] <- seq(compCond,NumCond*NumReps,NumCond)
              RR[2,seq(j,ncomps*NumReps,ncomps)] <- seq(refCond,NumCond*NumReps,NumCond)
            }
            Comps$RR <- RR

            compNames <- paste(allComps[,2],"vs",allComps[,1],sep="_")
            Comps$compNames <- compNames
            Comps$num <- nrow(allComps)
            Comps$allComps <- allComps
          }

          if(ncol(dat) == input$NumReps*input$NumCond) {
            updateCollapse(session,"Input",open = "Statistical testing")
          } else {
            updateCollapse(session,"Input",close = "Statistical testing")
          }

          paste(ifelse(mode(as.matrix(dat))!="numeric","<b>Wrong file format
                       /setup</b></br>",""),
                ifelse(ndatcol != NumReps*NumCond,"<b>Column number doesn't fit
                       with number of replicates and conditions!</b><br/>",""),
                "Number of features: ",nrow(dat),
                "<br/>Number of data columns in file:", ndatcol,
                "<br/>Percentage of missing values:",
                round(sum(is.na(dat))/nrow(dat)/ncol(dat)*100,digits = 2),"<br/>",
                paste("<i>Condition ",1:NumCond,":</i>", sapply(1:NumCond, function(x)
                  paste(colnames(dat)[(0:(NumReps-1))*NumCond+x],collapse=", ")),"<br/>",collapse="")
          )
        })

      })

      ## Put the data into a SummarizedExperiment object
      quantDataMatrix <- as.matrix(dat)
      # Create the SummarizedExperiment object
      fulldata <- SummarizedExperiment(assays = list(quant = quantDataMatrix))
      if (NumCond * NumReps == ncol(dat)) {
        sampleMetadata <- data.frame(Condition = rep(paste0("C", 1:NumCond), NumReps),
                                     Replicate = rep(1:NumReps, each=NumCond))
        fulldata <- SummarizedExperiment(assays = list(quant = quantDataMatrix),
                                         colData = sampleMetadata)
      }

      if (!is.null(addInfo))
        rowData(fulldata) <- data.frame(addInfo)
      # Adding metadata
      metadata(fulldata) <- list(
        NumReps = NumReps,
        NumCond = NumCond
      )

      currse(fulldata)
      currdata(dat)

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
            lapply(paste("C",rep(1:NumCond,NumReps),"Rep",rep(1:NumReps,each=NumCond),sep="_")[1:min(ncol(dat),NumCond*NumReps)],th,style="text-align: center;border-left:thin solid;")
          ),
          tr(
            th('Feature',style="text-align: center;"),
            if(!is.null(addInfo))
              lapply(colnames(addInfo),th,style="text-align: center;border-left:thin solid;"),
            lapply(colnames(dat)[1:min(ncol(dat),NumReps*NumCond)], th,style="text-align: center;border-left:thin solid;")
          )
        )
      ))
      outTable <- NULL
      if (!is.null(addInfo)) outTable <- data.frame(addInfo,dat[,1:min(ncol(dat),NumReps*NumCond)],stringsAsFactors = F)
      else outTable <- dat[,1:min(ncol(dat),NumReps*NumCond)]
      output$stat_table <- DT::renderDataTable(DT::datatable(outTable,
                                                             filter = list(position = 'top', clear = FALSE),colnames = c('model' = 1),
                                                             options = list(scrollX = TRUE,dom = 'Blfrtip',
                                                                            columnDefs = list(list(width = '20%', targets = colnames(data))),
                                                                            autoWidth=T),
                                                             class="compact",
                                                             container=sketch))

    }
  })

  observe({
    #dev.control(displaylist="enable")
    input$button
    print(input$button)
    isolate({
      se <- currse()
      dat <- NULL
      if (!is.null(se)) {
        dat <- assay(se)
      }
      if (input$button == currButton)
        return()
      if (!is.null(dat)) {
        print("Run tests")
        currButton <<- input$button
        NumReps <- metadata(se)$NumReps
        NumCond <- metadata(se)$NumCond
        isPaired <- input$is_paired
        allComps <- Comps$allComps

        RR <- Comps$RR
        compNames <- Comps$compNames

        if (ncol(dat) == NumReps*NumCond) {

          withProgress(message="Calculating ...", min=0,max=1, {

            #           output$messages <- renderText("running")
            incProgress(0.1, detail = paste("Running statistical tests"))
            dat [!is.finite(as.matrix(dat))] <- NA
            dat <<- dat
            MAData<-NULL
            # Rearranged dat (reference condition comes first)
            # UData <- NULL
            NumComps <- Comps$num
            if (NumComps > 0) {
              if (isPaired) {
                fulldata <- PolySTest_paired(se, allComps)
              } else {
                fulldata <- PolySTest_unpaired(se, allComps)
              }


              setProgress(0.8, detail = paste("Preparing data"))
              print("Preparing data")

              updateCheckboxGroupInput(session,"selTests",choices = testNames2, selected = "PolySTest")
              updateCheckboxGroupInput(session,"selComps",choices = compNames, selected=compNames[1])
              Qvalue <- as.matrix(rowData(fulldata)[, grep("^FDR", colnames(rowData(fulldata)))])
              Pvalue <- as.matrix(rowData(fulldata)[, grep("^p_values", colnames(rowData(fulldata)))])
              LogRatios <- as.matrix(rowData(fulldata)[, grep("^log_ratios", colnames(rowData(fulldata)))])
              colnames(Qvalue) <- paste("FDR",rep(testNames2,each=NumComps),rep(compNames,length(testNames2)), sep="_")
              FullReg(cbind(rowData(fulldata),assay(fulldata)))

              currse(fulldata)

              # Calculate best fc and qlim combination for all tests but the t-test
              setProgress(0.7, detail = paste("Calculating favorable FDR and fc thresholds"))
              tcomb <- FindFCandQlim(Qvalue[, grep("(^FDR_PolySTest_)|(^FDR_limma_)|(^FDR_Miss_Test_)",names(Qvalue))],
                                     LogRatios)
              setProgress(0.9, detail = paste("Creating figures"))


              # Set fold-change thresholds
              updateSliderInput(session,"fcval1",min=round(min(LogRatios,na.rm=T),1),
                                max=0, value=c(-tcomb[1]))
              updateSliderInput(session,"fcval2",min=0,
                                max=round(max(LogRatios,na.rm=T),1), value=c(tcomb[1]))

              # Set q-value threshold
              updateNumericInput(session,"qval",value=tcomb[2])

              output$table_stats <- renderText(paste("Number selected features:",length(input$stat_table_rows_selected)))

              # Set plot window sizes
              heightSize(max(200, 200 + 20 * length(input$stat_table_rows_selected)))
              session$sendCustomMessage(type = "pval_plot", message = Comps$num * 400)
              session$sendCustomMessage(type = "volc_plot", message = Comps$num * 400)


              # Arrange table header
              sketch = htmltools::withTags(table(
                class = 'display',
                thead(
                  tr(
                    th('',style="text-align: center;"),
                    if(!is.null(addInfo))
                      th(colspan = ncol(addInfo), 'Metadata',style="text-align: center;border-left:thin solid;"),
                    th(colspan = NumComps, 'log_ratios',style="text-align: center;border-left:thin solid;"),
                    th(colspan = NumComps*NumTests, 'FDRs',style="text-align: center;border-left:thin solid;")
                  ),
                  tr(
                    th('',style="text-align: center;"),
                    if(!is.null(addInfo))
                      th(colspan = ncol(addInfo), ''),
                    th(colspan = NumComps, '',style="text-align: center;border-left:thin solid;border-left:thin solid;"),
                    th(colspan = NumComps, 'PolySTest',style="text-align: center;border-left:thin solid;color: #AA3333;"),
                    th(colspan = NumComps, 'limma',style="text-align: center;border-left:thin solid;"),
                    th(colspan = NumComps, 'Miss test',style="text-align: center;border-left:thin solid;"),
                    th(colspan = NumComps, 'rank products',style="text-align: center;border-left:thin solid;"),
                    th(colspan = NumComps, 'permutation test',style="text-align: center;border-left:thin solid;"),
                    th(colspan = NumComps, 't-test',style="text-align: center;border-left:thin solid;")
                  ),
                  tr(
                    th('Feature',style="text-align: center;"),
                    if(!is.null(addInfo))
                      lapply(colnames(addInfo),th,style="text-align: center;border-left:thin solid;"),
                    lapply(rep(compNames,NumTests+1), th,style="text-align: center;border-left:thin solid;")
                  )
                )
              ))
              output$stat_table <- DT::renderDataTable({
                FullReg(as.data.frame(FullReg()))
                DT::datatable(as.data.frame(FullReg()),
                              filter = list(position = 'top', clear = FALSE),colnames = c('feature' = 1),
                              options = list(scrollX = TRUE,dom = 'Blfrtip',
                                             columnDefs = list(list(width = '20%', targets = colnames(FullReg()))),
                                             autoWidth=T,lengthMenu = c(5, 10, 50,100),
                                             buttons = list('colvis', 'copy', 'print')),
                              extensions=c("Buttons","FixedColumns"),class="compact",
                              container=sketch) %>%
                  formatSignif(grep("log_ratios",colnames(FullReg()),value=T),digits=2) %>%
                  formatSignif(grep("FDR",colnames(FullReg()),value=T),digits=2)

              })
              proxy <<- dataTableProxy('stat_table')

              observeEvent(input$resetSelection, {
                proxy %>% DT::selectRows(NULL)
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

              incProgress(0.9, detail = paste("Plotting first results"))

              print("Preparing data done")

              ## Plotting DRFs vs q-value thresholds + volcano plots + p-val vs q-vals
            }
          })
        }
      }
    })


  })

  ## Delay reaction to selecting rows in data table
  triggerUpdate <- debounce(reactive(list(input$stat_table_rows_selected,
                                          input$resetSelection,
                                          input$allSelection,
                                          input$allPageSelection,
                                          input$selComps,
                                          input$selTests)

  ),1000)


  observeEvent(input$allLimsSelection, {
    fulldata <- currse()
    if(!is.null(fulldata) & !is.null(FullReg())) {
      print("allLimsSelection started")
      selCols <- paste("FDR",apply(expand.grid(input$selTests, input$selComps), 1, paste, collapse="_"), sep="_")
      fclim <- c(input$fcval1, input$fcval2)
      NumComps <- Comps$num
      FCRegs <- filterFC(fulldata, NumTests, NumComps, fclim)
      selFeat <- which(rowSums(FCRegs[,which(colnames(FCRegs) %in% selCols),drop=F]<input$qval)>0)
      proxy %>% DT::selectRows(NULL)
      if (length(selFeat) > 0)
        proxy %>% DT::selectRows(as.numeric(selFeat))
      print("allLimsSelection finished")
    }
  })

  # volcano plots
  output$plotvolc <- renderPlot({
    triggerUpdate()
    if(!is.null(FullReg())) {
      fulldata <- currse()
      compNames <- Comps$compNames
      input$button
      input$stat_table
      print("plot volcano plots")
      sel_prots <- triggerUpdate()[[1]]
      qlim <- input$qval
      fclim <- c(input$fcval1, input$fcval2)
      plotVolcano(fulldata, compNames, testNames2, sel_prots, qlim, fclim)
      output$downloadVolcanoPdf <- downloadHandler(
        filename = function() {
          paste("VolcanoPlots", Sys.Date(), ".pdf", sep="");
        },
        content = function(file) {
          pdf(file,height=12,width=20)
          plotVolcano(fulldata, compNames, testNames2, sel_prots, qlim, fclim)
          dev.off()
        }
      )
    }
  })

  output$plotexpression <- renderPlot({
    triggerUpdate()
    fulldata <- currse()
    if(!is.null(fulldata) & !is.null(FullReg())) {
      compNames <- Comps$compNames
      print("plot expressions")
      input$button
      input$stat_table
      input$profiles_scale

      sel_prots <- triggerUpdate()[[1]]
      qlim <- input$qval
      fclim <- c(input$fcval1, input$fcval2)
      if (length(sel_prots)> 0) {
        withProgress(message="Creating expression profiles ...", min=0,max=1, {
          setProgress(0.5)

          # par(mfrow=c(1,3))
          plotExpression(fulldata, compNames, testNames2, sel_prots, input$profiles_scale, qlim, fclim)

          output$downloadExprPdf <- downloadHandler(
            filename = function() {
              paste("SelExpressionProfiles", Sys.Date(), ".pdf", sep="");
            },
            content = function(file) {
              pdf(file,height=8,width=25)
              plotExpression(fulldata, compNames, testNames2, selProts, input$profiles_scale, qlim, fclim)
              dev.off()

            })
        })
      }
    }
  },height=400)

  output$plotheatmap <- renderPlotly({
    triggerUpdate()
    fulldata <- currse()
    if(!is.null(fulldata) & !is.null(FullReg())) {
      compNames <- Comps$compNames
      NumComps <- Comps$num

      print("Plot heatmap")
      input$button
      input$stat_table
      input$heatmap_scale
      sel_prots <- triggerUpdate()[[1]]
      qlim <- input$qval
      fclim <- c(input$fcval, input$fcval2)
      isolate({
        p <- plotHeatmaply(fulldata, NumComps, sel_prots, input$heatmap_scale)
      })
      output$downloadHeatmapPdf <- downloadHandler(
        filename = function() {
          paste("Heatmap", Sys.Date(), ".pdf", sep="");
        },
        content = function(file) {
          plotHeatmaply(fulldata, NumComps, sel_prots, input$heatmap_scale, file=file)
        })
      p
    }
  })#,height=800)

  output$plotregdistr <- renderPlot({
    triggerUpdate()
    fulldata <- currse()
    if(!is.null(fulldata) & !is.null(FullReg())) {
      compNames <- Comps$compNames
      NumComps <- length(compNames)
      NumTests <- length(testNames2)

      print("Plot q-value numbers")
      qlim <- input$qval
      fclim <- c(input$fcval1, input$fcval2)
      input$button
      isolate({
        print(plotUpset(fulldata, NumTests, NumComps, qlim, fclim))
        output$downloadUpSetPdf <- downloadHandler(
          filename = function() {
            paste("UpSetProfiles", Sys.Date(), ".pdf", sep="");
          },
          content = function(file) {
            pdf(file,height=8,width=8)
            print(plotUpset(fulldata, NumTests, NumComps, qlim, fclim))
            dev.off()
          })
      })
    }

  }, height=400)

  output$plotreg <- renderPlot({
    fulldata <- currse()
    triggerUpdate()

    if(!is.null(fulldata) & !is.null(FullReg())) {
      NumComps <- length(Comps$compNames)

      print("Plot q-value number")
      qlim <- input$qval
      fclim <- c(input$fcval1, input$fcval2)
      input$button
      plotRegNumber(fulldata, NumComps, testNames2, qlim, fclim)
      output$downloadRegDistrPdf <- downloadHandler(
        filename = function() {
          paste("RegDistrPlots", Sys.Date(), ".pdf", sep="");
        },
        content = function(file) {
          pdf(file,height=8,width=8)
          plotRegNumber(fulldata, NumComps, testNames2, qlim, fclim)
          dev.off()
        })
    }

  },height=400)

  output$downloadData <- downloadHandler(
    filename = function() {
      paste("Results", Sys.Date(), ".csv", sep="");
    },
    content = function(file) {
      fulldata <- currse()
      if(!is.null(fulldata) & !is.null(FullReg())) {
        write.csv(cbind(rowData(fulldata), assay(fulldata)), file)
      } else {
        showNotification(
          "Warning: No results yet. Please ensure all required fields are filled and
          the analysis is done.",
          type = "warning", # Notification type: "default", "message", "warning", or "error"
          duration = 5 # How long to show the notification, in seconds
        )
      }
    })

  # Send data back to calling OmicsQ
  observeEvent(input$retrieve_output, isolate({
    #print(input$retrieve_output)
    # make table in right format
    if (input$retrieve_output == "Get data" & !is.null(fulldata) & !is.null(FullReg())) {
      print("Sending data back")
      fulldata <- currse()
      outdata <- as.data.frame(rowData(fulldata))
      print(head(outdata, 1))
      BackMessage <- toJSON(list(expr_matrix=as.list(outdata)))
      js$send_results(dat=BackMessage)
    }
  }))

  output$plotpval <- renderPlot({
    proxy
    if(!is.null(FullReg())) {
      print("plotting pvalue distributions")
      fulldata <- currse()
      compNames <- Comps$compNames
      output$downloadPvalueDistrPdf <- downloadHandler(
        filename = function() {
          paste("PvalueDistrPlots", Sys.Date(), ".pdf", sep="");
        },
        content = function(file) {
          pdf(file,height=12,width=20)
          plotPvalueDistr(fulldata, compNames, testNames)
          dev.off()
        }
      )
      plotPvalueDistr(fulldata, compNames, testNames)
    }
  })


})
