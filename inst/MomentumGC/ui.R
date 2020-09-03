# This is a Shiny web application. You can run the application by executing `vngGrowthCurves:::Momentum_GC_run()`

library(shiny)

# Define UI for application that plot growth curves
shinyUI(
  fluidPage(
  
  plotly::plotlyOutput('gc_plot', height="auto"),
  
  hr(),
  
  h3("Momentum growth curves"),
  fluidRow(
    tags$head(
      tags$style(type="text/css", "label{ display: table-cell; text-align: center; vertical-align: middle; padding-right: 0.5em;}
            .form-group { display: table-row;}")
    ),

    column(7, # files selection
           # fluidRow(
           #   column(3, shinyFiles::shinyDirButton("path", "Select directory", "Select", buttonType="primary")),
           #   column(9, verbatimTextOutput("path")),
           # ),
           splitLayout(
             cellWidths = c("9em", NA),
             shinyFiles::shinyDirButton("path", "Select directory", "Select", buttonType="primary"),
             verbatimTextOutput("path")
           ),
           fluidRow(
             column(9, textInput("re", "Filter", 
                                 placeholder = "Path pattern (regexp)")),
             column(2, actionButton("refreshPlot","Plot ", icon=icon("chart-line"), class = "btn-primary"),
                    # style="float:right;"),
                    # class="pull-sm-right"), # will work the day shiny uses BS4
                    ),
           )
    ),
    
    column(4, offset = 1, # options
                  textInput("range", "Wells range", 
                     placeholder = "e.g. A1:D8 (empty for all wells)"),
           
                  checkboxInput('log', 'log scale'),
                  numericInput('blank', 'Blank value', value=0, min=0, max=1, step=.001),
    )
  ),
  
  fluidRow(
    
    # column(7, verbatimTextOutput("series", placeholder=FALSE) ),
    column(7, wellPanel(htmlOutput("series"))),
    # column(7, conditionalPanel("output.path_set == 'yes'", wellPanel(htmlOutput("series")))),
  )
  )
)