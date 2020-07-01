library(shinyBS)
# library(d3heatmap)
library(heatmaply)
library(shinydashboard)
title <- tags$img(src="Logo.svg",style="width:200px")

shinyUI(dashboardPage(skin="blue",
                      dashboardHeader(title=title,dropdownMenu(    
                        type = "notifications", 
                        icon = icon("question-circle"),
                        badgeStatus = NULL,
                        headerText = "See also:",
                        
                        notificationItem("source code and installation instructions", icon = icon("file"),
                                         href = "http://bitbucket.com/veitveit/polystest"),
                        notificationItem("paper", icon = icon("scroll"),
                                         href = "https://doi.org/10.1074/mcp.RA119.001777"),
                        notificationItem("author", icon = icon("address-card"),
                                         href = "http://computproteomics.bmb.sdu.dk"),
                        notificationItem("institution", icon = icon("university"),
                                         href = "sdu.dk")
                      )
                      ),
                      dashboardSidebar(width = "300",
                                       tags$style(HTML("
      .panel-primary {
      background-color: #333;
      }
      section.sidebar .shiny-input-container {
           padding: 0px;
      }

    ")
                                       ),
                                       fluidPage(
                                       h5("A tool to determine and visualize differentially regulated features using multiple approaches. For more information, press the question mark button on the upper right."),
                                       bsCollapse(id="Input",multiple=T,open="Data input",
                                                  bsCollapsePanel("Data input",style="primary",
                                                                  fileInput("in_file", "Input file:",accept=c("txt/csv", "text/comma-separated-values,text/plain", ".csv")),
                                                                  actionLink("example","Load example"),
                                                                  bsTooltip("example","Example data set from proteomics data (4 conditions, 3 replicates)",trigger="hover"),
                                                                  fluidRow(
                                                                    checkboxInput(inputId="is_header", label="Column names?", value=TRUE),
                                                                    bsTooltip("is_header","Input file has column headers (only first row)",trigger="hover"),
                                                                    checkboxInput(inputId="row.names", label="Row names?", value=TRUE),
                                                                    bsTooltip("row.names","First column contains unique feature names",trigger="hover"),
                                                                    radioButtons(inputId="delimiter", "Delimiter character", choices=c(",",";","tab"), selected=",",inline=T),
                                                                    bsTooltip("delimiter","Cells are separated by ...",trigger="hover"),
                                                                    radioButtons(inputId="digits", "Decimal character", choices=c(".",","), selected=".",inline=T),
                                                                    bsTooltip("digits","Used character for digits",trigger="hover")
                                                                  )
                                                  ),
                                                  bsCollapsePanel("Data layout",style="primary",
                                                                  numericInput("ColQuant",min=2,max=20,value=2,label="First column for quantification",step=1),
                                                                  bsTooltip("ColQuant","Number of first column that contains the to-be-analyzed values (e.g. 3 when the first two column contain protein IDs and protein descriptions, respectively)",trigger="hover"),
                                                                  checkboxInput(inputId="qcol_order", label="Replicates are grouped",value = T),
                                                                  bsTooltip("qcol_order","Given that replicates are numbered and conditions as given by letters A,B,... , grouped replicates denotes A1,B1,C1,...,A2,B2,C2 ..., while ungrouped means A1,A2,A3,...,B1,B2,B3,...",trigger="hover"),
                                                                  numericInput("NumReps",min=2,max=20,value=3,label="Number of replicates",step=1),
                                                                  bsTooltip("NumReps","Number of replicates per condition (fill by empty columns when different for different conditions)",trigger="hover"),
                                                                  numericInput("NumCond",min=2,max=20,value=2,label="Number of conditions",step=1),
                                                                  bsTooltip("NumCond","Number of experimental conditions that should be compared",trigger="hover")
                                                  ),
                                                  # numericInput("refCond",min=1,max=20,value=1,label="Reference condition",step=1),
                                                  # bsTooltip("refCond","Experimental condition to which the other conditions will be compared to (e.g. 1 vs 2, 3 vs 2, 4 vs 2 - when set to 2)",trigger="hover"),
                                                  bsCollapsePanel("Statistical testing",style="primary",
                                                                  actionButton("button","Run analysis",icon = icon("play"),style="background-color:#FFED66"),
                                                                  bsTooltip("button","Run statistical tests and their evaluation",trigger="hover"),
                                                                  checkboxInput(inputId="is_paired", label="Paired tests?", value=F),
                                                                  bsTooltip("is_paired","Paired: tests are paired within the replicates; Unpaired: comparison each versus each",trigger="hover"),
                                                                  uiOutput("stat_comparisons"),
                                                                  actionButton("addComp", "Add new comparison"),
                                                                  bsTooltip("addComp","Add new pair of conditions to be tested for differentially regulated features")
                                                  )))),
                      dashboardBody(
                        h3(""),
                        box(title="Help (click on the right to see the help page)",width=12, solidHeader = T, status="warning",collapsible = T,collapsed=T,
                            htmlOutput("description"),
                            htmlOutput("description2"),
                            htmlOutput("description3"),
                            htmlOutput("description4")),
                        br(),
                        box(title="Data details", width=3,solidHeader = T,status="info",
                            htmlOutput("input_stats")),
                        conditionalPanel(
                          condition = "input.button > 0",
                          box(title="Thresholds",width=3,solidHeader = T,status="primary",
                              p("To select feature groups and adapt visualization groups below."),
                              numericInput("qval","FDR threshold",value=0.01,min=0,max=1,step=0.01,width = "200px"),
                              bsTooltip("qval","FDR threshold to assign significant features",trigger="hover"),
                              fluidRow(
                                column(6,
                              numericInput("fcval1","lower log-ratio threshold", value=0, min=-1, max=0, step=0.1 )),
                              column(6,
                              numericInput("fcval2","upper log-ratio threshold", value=0, min=0, max=1, step=0.1 ))),
                              bsTooltip("fcval1","Filter table for log-ratio (log of fold change) threshold",trigger="hover"),
                              bsTooltip("fcval2","Filter table for log-ratio (log of fold change) threshold",trigger="hover"),
                              actionButton("allLimsSelection",HTML("Select all features filtered by<br/> FDR and fold-change"))),
                          box(title="Select tests and comparison",width=2,solidHeader = T,status="primary",
                              p("Further criteria for feature selection."),
                              checkboxGroupInput("selTests",label = "Which statistical test(s)?",choices = list("Please wait ...")),
                              bsTooltip("selTests", "Take only features that fulfill the conditions above in the given statistical tests", trigger="hover"),
                              checkboxGroupInput("selComps",label = "Which comparison(s)?",choices = list("Please wait ...")),
                              bsTooltip("selComps", "Take only features that fulfill found to be regulated in the following comparisons", trigger="hover")),
                          box(title="Select features",width=3,solidHeader = T,status="primary",
                              htmlOutput("table_stats"),
                              actionButton("allPageSelection","Select all shown features"),br(),
                              actionButton("allSelection","Select all (filtered) features"),br(),
                              actionButton("resetSelection","Clear selected features in table"),br(),br(),
                              downloadButton('downloadData', 'Download results',style="background-color:#FFAAAA"),
                              bsTooltip("downloadData","Download the entire data table with log-ratios and FDRs",trigger="hover"))
                          ,br(),hr()),
                        # conditionalPanel("$('#dtable_out').hasClass('recalculating')",tags$div('Loading ... ')),
                        box(status="info",solidHeader=T,width=12,
                            column(div(DT::dataTableOutput("stat_table"),style="font-size:100%"),width=12)),
                        conditionalPanel(
                          condition = "input.button > 0",
                          box(title="Details on tests and feature (e.g. proteins) expression changes (max. 30 with lowest unified FDRs shown)",
                              collapsible = TRUE,status="success",solidHeader = T,collapsed=T,
                              checkboxInput("profiles_scale","Scale features to mean"),
                              plotOutput("plotexpression",height="auto"),
                              downloadButton("downloadExprPdf","Download as pdf"),width=12),
                          box(title="Co-expression patterns and significance",collapsible = TRUE,status="success",solidHeader = T,collapsed=T,
                              # d3heatmapOutput("plotheatmap",height="auto"),width=3),
                              # plotOutput("plotheatmap",height="auto")
                              checkboxInput("heatmap_scale","Scale features to mean"),
                              plotlyOutput("plotheatmap", height = "700px"),
                              downloadButton("downloadHeatmapPdf","Download as pdf")
                              ,width=12),
                          box(title="Comparison of tests and conditions (volcano plots)",collapsible = TRUE,status="success",solidHeader = T,
                              plotOutput("plotvolc",height="auto"),
                              downloadButton("downloadVolcanoPdf","Download as pdf")
                              ,width=12),
                          box(title="Feature numbers per test and condition",
                              collapsible = TRUE,status="success",solidHeader = T,
                              p("This figure compares the number of regulated features for all statistical tests across comparisons. 
                                This applies the given thresholds for FDR and log-ratios above. Here, you can compare overlap between
                                performance and between conditions."),
                              plotOutput("plotregdistr",click="plotregdistr_click",height="auto"),
                              downloadButton("downloadUpSetPdf","Download as pdf")
                              ,width=12),
                          
                          box(title="Number of significant features versus thresholds",collapsible = TRUE,status="success",solidHeader = T,
                              plotOutput("plotreg",height="auto"),
                              downloadButton("downloadRegDistrPdf","Download as pdf")
                              ,width=12)),
                        box(title="Distribution of p-values (not corrected for multiple testing)",collapsible = TRUE,status="success",solidHeader = T,
                            plotOutput("plotpval",height="auto"),
                            downloadButton("downloadPvalueDistrPdf","Download as pdf")
                            ,width=12)
                        
                      )
)
)


