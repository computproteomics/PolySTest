library(shinyBS)
library(shinydashboard)
shinyUI(dashboardPage(
  dashboardHeader(title="PolySTest"),
  dashboardSidebar(width = "300",
                   fluidPage(
                     h5("A tool to determine differentially regulated features using multiple approaches"),
                     fileInput("in_file", "Input file:",accept=c("txt/csv", "text/comma-separated-values,text/plain", ".csv")),
                     actionButton("button","Run analysis"),
                     actionLink("example","Load example"),
                     bsTooltip("example","Example data set from proteomics data (4 conditions, 3 replicates)",trigger="hover"),
                     checkboxInput(inputId="is_paired", label="Paired tests?", value=T),
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
    htmlOutput("description"),
    htmlOutput("description2"),
    htmlOutput("description3"),
    htmlOutput("description4"),
    br(),
      fluidRow(
        
        box(title="Data details", 
            fluidRow(
              column(width=6,htmlOutput("input_stats")),
              conditionalPanel(
                condition = "input.button > 0",
                column(width=6,numericInput("qval","q-value threshold",value=0.01,min=0,max=1),
                     bsTooltip("qval","q-value threshold to assign significant features",trigger="hover")
              ),
              column(width=6,
                     actionButton("resetSelection","Clear selected rows in table"),
                     downloadButton('downloadData', 'Download results',style="background-color:#FFAAAA"),
                     bsTooltip("downloadData","Download the entire data table with log-ratios and q-values",trigger="hover")
              ))),br(),hr(),
            # conditionalPanel("$('#dtable_out').hasClass('recalculating')",tags$div('Loading ... ')),
            column(DT::dataTableOutput("stat_table"),width=12),width = 12),
        box(title="Volcano plots",
            plotOutput("plotvolc",height="auto"),width=12),
        box(title="p-value histograms",
            plotOutput("plotpval",height="auto")),
        box(title="Expression profiles",
            plotOutput("plotexpression",height="auto")),
        box(title="Distribution regulated features over different tests and conditions",
            plotOutput("plotregdistr",click="plotregdistr_click",height="auto"),width=12),
        box(title="Percentage of regulated features",
            plotOutput("plotreg",height="auto"),width=12)
      )
    )
  
)
)
