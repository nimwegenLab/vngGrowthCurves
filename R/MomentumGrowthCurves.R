#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/

globalVariables(c("data", "datetime", "filedate", "filename", "filetime", "is_txt", "path", "series", "time_sec"))


range_to_wells <- function(c_lim) {
  # NB: excel cells use letters for columns while microplate wells use them for rows... 
  tidyr::expand_grid(col=c_lim$ul[1]:c_lim$lr[1], row=LETTERS[c_lim$ul[2]:c_lim$lr[2]])
}

# Define UI for application that plot growth curves
Momentum_GC_ui <- shiny::fluidPage(
  
  plotly::plotlyOutput('gc_plot', height="720px"),
  
  shiny::hr(),
  
  shiny::h3("Momentum growth curves"),
  shiny::fluidRow(
    shiny::tags$head(
      shiny::tags$style(type="text/css", "label{ display: table-cell; text-align: center; vertical-align: middle; padding-right: 0.5em;}
            .form-group { display: table-row;}")
    ),

    shiny::column(7, # files selection
           # fluidRow(
           #   column(3, shinyFiles::shinyDirButton("path", "Select directory", "Select", buttonType="primary")),
           #   column(9, verbatimTextOutput("path")),
           # ),
           shiny::splitLayout(
             cellWidths = c("9em", NA),
             shinyFiles::shinyDirButton("path", "Select directory", "Select", buttonType="primary"),
             shiny::verbatimTextOutput("path")
           ),
           shiny::fluidRow(
             shiny::column(9, shiny::textInput("glob", "Pattern", 
                                 placeholder = "Path pattern (globbing accepted)")),
             shiny::column(2, shiny::actionButton("refreshPlot","Plot ", icon=shiny::icon("chart-line"), class = "btn-primary"),
                    # style="float:right;"),
                    # class="pull-sm-right"), # will work the day shiny uses BS4
                    ),
           )
    ),
    
    shiny::column(4, offset = 1, # options
                  shiny::textInput("range", "Wells range", 
                     placeholder = "e.g. A1:D8 (empty for all wells)"),
           
                  shiny::checkboxInput('log', 'log scale'),
                  shiny::numericInput('blank', 'Blank value', value=0, min=0, max=1, step=.001),
    )
  ),
  
  shiny::fluidRow(
    
    # column(7, verbatimTextOutput("series", placeholder=FALSE) ),
    shiny::column(7, shiny::wellPanel(shiny::htmlOutput("series"))),
    # column(7, conditionalPanel("output.path_set == 'yes'", wellPanel(htmlOutput("series")))),
  )
  )

# Define server logic required to draw a histogram
Momentum_GC_server <- function(input, output) {

  # enable directory selection in UI
  volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), shinyFiles::getVolumes()())
  shinyFiles::shinyDirChoose(
    input,
    'path',
    roots = volumes,
  )
  output$path <- shiny::renderPrint({
    if (is.integer(input$path)) {
      cat("No input directory has been selected.")
    } else {
      shinyFiles::parseDirPath(volumes, input$path)
    }
  })
  
  myseries <- shiny::reactive({
    if (is.integer(input$path)) {
      NULL
    } else {
      glob <- if (input$glob=="") NULL else input$glob
      # browser()
      shinyFiles::parseDirPath(volumes, input$path) %>% 
        fs::dir_ls(type="file", glob=glob) %>%
        dplyr::tibble(path=.) %>%
        dplyr::mutate(filename = fs::path_file(path),
               is_txt=stringr::str_detect(filename, "\\.txt$"),
               series = stringr::str_extract(filename, "^.*_pl\\d+")
        ) %>%
        dplyr::group_by(series) %>%
        tidyr::nest() %>%
        identity()
    }
  })
  
  output$series <- shiny::renderPrint({
    if (is.null(myseries())) {
      base::invisible() # return nothing (as per the doc)
    } else {
      if (any(is.na(myseries()$series)))
        print("Message: found files not matching the pattern...")
      
      if (myseries() %>% dplyr::filter(is.na(series)) %>% nrow())
        if (myseries() %>% dplyr::filter(is.na(series)) %>% tidyr::unnest(data) %>% dplyr::filter(!is_txt) %>% nrow())
          print("Warning: ignoring files with extension different from `.txt`...")
      
      cat("Found ", myseries() %>% dplyr::filter(!is.na(series)) %>% nrow(), " files series:")
      if (myseries() %>% dplyr::filter(!is.na(series)) %>% nrow())
        myseries() %>% dplyr::filter(!is.na(series)) %>% 
        tidyr::unnest(data) %>% dplyr::filter(is_txt) %>% 
        dplyr::group_by(series) %>% dplyr::summarise(n_files = dplyr::n(), .groups="drop_last") %>%
        knitr::kable(format='html')
    }
  })

  parse_data <- shiny::reactive({
    if (!is.null(myseries()))
      if (myseries() %>% dplyr::filter(!is.na(series)) %>% nrow()) {
        mydata <- myseries() %>% dplyr::filter(!is.na(series)) %>% 
          tidyr::unnest(data) %>% dplyr::filter(is_txt) %>% 
          tidyr::extract(filename, c('filedate', 'filetime'), "_pl1_(\\d{8})_(\\d{6})(?:_FE_)?.txt$", remove=FALSE) %>% 
          dplyr::rowwise() %>% 
          dplyr::mutate(
            datetime=lubridate::ymd_hms(paste0(filedate, filetime)), filedate=NULL, filetime=NULL,
            data=list(vngGrowthCurves::read_Biotek_Synergy2_matrix(path)),
            # data=purrr::map(path, ~vngGrowthCurves::read_Biotek_Synergy2_matrix(.))
          ) %>%
          tidyr::unnest(data)
      }
  })
        
  create_plot <- shiny::eventReactive(c(input$refreshPlot, input$log), {
    if (!is.null(myseries()))
      if (myseries() %>% dplyr::filter(!is.na(series)) %>% nrow()) {
        mydata <- parse_data()
        
        c_lim <- try(cellranger::as.cell_limits(input$range), silent = TRUE)
        if (class(c_lim)[1] != 'try-error')
          mydata <- dplyr::semi_join(mydata, range_to_wells(c_lim), by=c("row", "col"))
        
        if (input$blank > 0)
          mydata <- dplyr::mutate(mydata, value=value-input$blank)
          
        pl <- mydata %>% 
          dplyr::group_by(series, row, col) %>%
          dplyr::mutate(time_sec=as.numeric(datetime - min(datetime))) %>%
          # identity() %>% print()
          ggplot2::ggplot() +
          ggplot2::facet_grid(row~col) +
          ggplot2::geom_point(ggplot2::aes(time_sec / 3600, value, col=series), size=1, stroke=0, alpha=.7) +
          ggplot2::labs(x="time (h)") +
          ggplot2::theme_bw() +
          ggplot2::theme(legend.position = 'bottom') +
          NULL
        
        if (input$log)
          pl <- pl + ggplot2::scale_y_log10()
        
        return(pl)
      }
  })
  
  output$gc_plot <- plotly::renderPlotly({
    if (!is.null(create_plot()))
        plotly::ggplotly(create_plot()) %>%
          plotly::layout(legend = list(orientation = "h", x=0, y=-0.07, itemsizing="constant"))
  })      
  
}

# Run the application 
Momentum_GC_run <- function()
  shiny::shinyApp(Momentum_GC_ui, Momentum_GC_server)
