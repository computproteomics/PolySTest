library(shinyBS)
# library(d3heatmap)
library(heatmaply)
library(shinydashboard)
shinyUI(dashboardPage(skin="blue",
                      dashboardHeader(title="PolySTest"),
                      dashboardSidebar(width = "300",
                                       fluidPage(
                                         h5("A tool to determine differentially regulated features using multiple approaches"),
                                         fileInput("in_file", "Input file:",accept=c("txt/csv", "text/comma-separated-values,text/plain", ".csv")),
                                         actionButton("button","Run analysis"),
                                         bsTooltip("button","Run statistical tests and their evaluation",trigger="hover"),
                                         actionLink("example","Load example"),
                                         bsTooltip("example","Example data set from proteomics data (4 conditions, 3 replicates)",trigger="hover"),
                                         checkboxInput(inputId="is_paired", label="Paired tests?", value=F),
                                         bsTooltip("is_paired","Paired: tests are paired within the replicates; Unpaired: comparison each versus each",trigger="hover"),
                                         h3("Input file parameters"),
                                         fluidRow(
                                           checkboxInput(inputId="is_header", label="Column names?", value=TRUE),
                                           bsTooltip("is_header","Input file has column headers (only first row)",trigger="hover"),
                                           checkboxInput(inputId="row.names", label="Row names?", value=TRUE),
                                           bsTooltip("row.names","First column contains unique feature names",trigger="hover"),
                                           radioButtons(inputId="delimiter", "Delimiter character", choices=c(",",";","tab"), selected=",",inline=T),
                                           bsTooltip("delimiter","Cells are separated by ...",trigger="hover"),
                                           radioButtons(inputId="digits", "Decimal character", choices=c(".",","), selected=".",inline=T),
                                           bsTooltip("digits","Used character for digits",trigger="hover")
                                         ),
                                         h3("Data layout"),
                                         sliderInput("ColQuant",min=1,max=20,value=2,label="First column for quantification",step=1),
                                         bsTooltip("ColQuant","Number of first column that contains the to-be-analyzed values (e.g. 3 when the first two column contain protein IDs and protein descriptions, respectively)",trigger="hover"),
                                         sliderInput("NumReps",min=2,max=20,value=3,label="Number of replicates",step=1),
                                         bsTooltip("NumReps","Number of replicates per condition (fill by empty columns when different for different conditions)",trigger="hover"),
                                         sliderInput("NumCond",min=2,max=20,value=4,label="Number of conditions",step=1),
                                         bsTooltip("NumCond","Number of experimental conditions that should be compared",trigger="hover"),
                                         sliderInput("refCond",min=1,max=20,value=1,label="Reference condition",step=1),
                                         bsTooltip("refCond","Experimental condition to which the other conditions will be compared to (e.g. 1 vs 2, 3 vs 2, 4 vs 2 - when set to 2)",trigger="hover")
                                       )),
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
                              numericInput("qval","q-value threshold",value=0.01,min=0,max=1,step=0.01,width = "200px"),
                              bsTooltip("qval","q-value threshold to assign significant features",trigger="hover"),
                              sliderInput("fcval","fold-change threshold",value=c(-1,1),min=-1,max=1,step=0.1,width = "200px"),
                              bsTooltip("fcval","Filter table for fold-change threshold",trigger="hover"),
                              actionButton("allLimsSelection",HTML("Select all features filterd by<br/> q-value and fold-change"))),
                          box(title="Select tests and comparison",width=2,solidHeader = T,status="primary",
                              checkboxGroupInput("selTests",label = "Which statistical test(s)?",choices = list("Initiate")),
                              bsTooltip("selTests", "Take only features that fulfill the conditions above in the given statistical tests", trigger="hover"),
                              checkboxGroupInput("selComps",label = "Which comparison(s)?",choices = list("Initiate")),
                              bsTooltip("selComps", "Take only features that fulfill found to be regulated in the following comparisons", trigger="hover")),
                          box(title="Select features",width=3,solidHeader = T,status="primary",
                              htmlOutput("table_stats"),
                              actionButton("allPageSelection","Select all shown features"),br(),
                              actionButton("allSelection","Select all (filtered) features"),br(),
                              actionButton("resetSelection","Clear selected features in table"),br(),br(),
                              downloadButton('downloadData', 'Download results',style="background-color:#FFAAAA"),
                              bsTooltip("downloadData","Download the entire data table with log-ratios and q-values",trigger="hover"))
                          ,br(),hr()),
                          # conditionalPanel("$('#dtable_out').hasClass('recalculating')",tags$div('Loading ... ')),
                          box(status="info",solidHeader=T,width=12,
                              column(div(DT::dataTableOutput("stat_table"),style="font-size:100%"),width=12)),
                        conditionalPanel(
                          condition = "input.button > 0",
                          box(title="Expression profiles (max. 30 with lowest unified q-values shown)",
                              collapsible = TRUE,status="success",solidHeader = T,collapsed=T,
                              plotOutput("plotexpression",height="auto"),
                              downloadButton("downloadExprPdf","Download as pdf"),width=12),
                          box(title="Clustered data",collapsible = TRUE,status="success",solidHeader = T,collapsed=T,
                              # d3heatmapOutput("plotheatmap",height="auto"),width=3),
                              # plotOutput("plotheatmap",height="auto")
                              plotlyOutput("plotheatmap", height = "700px")
                              ,width=12),
                          box(title="Volcano plots",collapsible = TRUE,status="success",solidHeader = T,
                              plotOutput("plotvolc",height="auto"),
                              downloadButton("downloadVolcanoPdf","Download as pdf")
                          ,width=12),
                          box(title="Distribution regulated features over different tests and conditions",
                              collapsible = TRUE,status="success",solidHeader = T,
                              plotOutput("plotregdistr",click="plotregdistr_click",height="auto"),
                              downloadButton("downloadUpSetPdf","Download as pdf")
                          ,width=12),
                          
                          box(title="Percentage of regulated features",collapsible = TRUE,status="success",solidHeader = T,
                              plotOutput("plotreg",height="auto"),
                              downloadButton("downloadRegDistrPdf","Download as pdf")
                              ,width=12)),
                          box(title="p-value histograms",collapsible = TRUE,status="success",solidHeader = T,
                              plotOutput("plotpval",height="auto"),
                              downloadButton("downloadPvalueDistrPdf","Download as pdf")
                              ,width=12)
                        
                      )
)
)


