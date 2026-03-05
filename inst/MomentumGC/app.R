# This is a Shiny web application. You can run the application by executing `vngGrowthCurves:::Momentum_GC_run()`

library(shiny)

# library(shinyFiles)
library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)
library(vngGrowthCurves)

# =========================
# UI
# =========================

ui <- fluidPage(
  
  plotly::plotlyOutput("gc_plot", height="auto"),
  
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
             column(2,
                    actionButton("refreshPlot","Plot ", icon=icon("chart-line"), class="btn-primary"),
                    # style="float:right;"),
                    # class="pull-sm-right"), # will work the day shiny uses BS4
             )
           )
    ),
    
    column(4, offset = 1, # options
           
           selectInput("channel", "Channel", 'No choices here yet'),
           
           textInput("range", "Wells range",
                     placeholder = "e.g. A1:D8 (empty for all wells)"),
           
           checkboxInput('log', 'log scale'),
           
           numericInput('blank', 'Blank value', value=0, min=0, max=1, step=.001)
    )
  ),
  
  fluidRow(
    # column(7, verbatimTextOutput("series", placeholder=FALSE) ),
    column(7, wellPanel(htmlOutput("series"))),
    # column(7, conditionalPanel("output.path_set == 'yes'", wellPanel(htmlOutput("series")))),
  )
)

# =========================
# SERVER
# =========================

server <- function(input, output, session) {
  
  # IMPORTANT! close session when running standalone
  # this is needed to terminate the R process when the
  # shiny app session ends. Otherwise, you end up with a zombie process
  if (!interactive()) {
    session$onSessionEnded(function() {
      stopApp()
      q("no")
    })
  }
  
  # DIRECTORY SELECTION ####  
  volumes <- c(
    Home = fs::path_home(),
    # "R Installation" = R.home(),
    shinyFiles::getVolumes()()
  )
  dir <- getOption("MomentumGCdir")
  if (!is.null(dir)) volumes <- c(rlang::set_names(fs::path_tidy(dir), fs::path_file(dir)), volumes)
  
  shinyFiles::shinyDirChoose(
    input,
    "path",
    roots = volumes,
    session = session
  )
  
  # chatGPT trick 1: isolate directory parsing to avoid unnecessary reactivity
  selected_dir <- reactive({
    req(!is.integer(input$path))
    isolate(shinyFiles::parseDirPath(volumes, input$path))
  })
  
  output$path <- renderPrint({
    if (is.integer(input$path)) {
      cat("No input directory has been selected.")
    } else {
      selected_dir()
    }
  })
  
  # SERIES DISCOVERY ####
  myseries <- reactive({
  # myseries is NULL when no file can be parsed
    if (is.integer(input$path)) return(NULL)
    
    re <- if (input$re=="") NULL else input$re
    .paths <- tryCatch(
      fs::dir_ls(selected_dir(), type="file", regexp=re), #dir_ls() crashes when regex is not correct
      error=function(e) character(0)
    )
    
    tibble(path=.paths) %>%
      mutate(
        filename = fs::path_file(path),
        is_txt = stringr::str_detect(filename, "\\.txt$"),
        series = stringr::str_extract(filename, "^.*_pl\\d+")
      ) %>%
      group_by(series) %>%
      nest() %>%
      mutate(
        channels = map(data, ~try(read_Biotek_Synergy2_matrices(.$path[.$is_txt][1], .ch_only=T)))
      )
  })
  
  output$series <- renderPrint({
    series <- myseries()
    
    if (is.null(series)) {
      invisible()
    } else if (series %>% pull(channels) %>% unique() %>% length() > 1) {
      print("Error: found files with different measurement channels, please adjust your selection criteria.")
    } else {     
      if (series %>% filter(is.na(series)) %>% nrow())
        print("Message: found files not matching the pattern...")

      if (series %>% filter(is.na(series)) %>% nrow())
        if (series %>% filter(is.na(series)) %>% unnest(data) %>% filter(!is_txt) %>% nrow())
          print("Warning: ignoring files with extension different from `.txt`...")
      
      cat("Found ", series %>% filter(!is.na(series)) %>% nrow(), " files series:")
      
      if (series %>% filter(!is.na(series)) %>% nrow())
        series %>%
        filter(!is.na(series)) %>%
        unnest(data) %>%
        filter(is_txt) %>%
        group_by(series) %>%
        summarise(n_files = n(), .groups="drop_last") %>%
        knitr::kable(format='html')
    }
  })
  
  # HELPERS ####
    
  range_to_wells <- function(c_lim) {
    # NB: excel cells use letters for columns while microplate wells use them for rows...
    expand_grid(col=c_lim$ul[1]:c_lim$lr[1], row=LETTERS[c_lim$ul[2]:c_lim$lr[2]])
  }
  
  parse_data <- eventReactive(myseries(), {
    series <- myseries()
    req(series)
    
    if (series %>% filter(!is.na(series)) %>% nrow()) {
      series %>%
        filter(!is.na(series)) %>%
        select(-channels) %>%
        unnest(data) %>%
        filter(is_txt) %>%
        extract(filename, c('filedate','filetime'),
                "_pl\\d+_(\\d{8})_(\\d{6})(?:_FE_)?.txt$",
                remove=FALSE) %>%
        mutate(
          datetime=lubridate::ymd_hms(paste0(filedate,filetime)),
          filedate=NULL,
          filetime=NULL
        ) %>%
        mutate(
          data = purrr::map(path, vngGrowthCurves::read_Biotek_Synergy2_matrices)
        ) %>%
        unnest(data) %>%
        unnest(data) %>%
        mutate(value=as.numeric(value))
    }
  })
  
  observeEvent(myseries(), {
    series <- myseries()
    req(series)
    
    .ch <- series %>% pull(channels) %>% unique()
    
    # chatGPT trick 2: update selector safely
    if (length(.ch)==1 && is.character(.ch[[1]]) && length(.ch[[1]])>0) {
      updateSelectInput(session, "channel",
        label="Channel",
        choices=.ch[[1]],
        selected=.ch[[1]][1]
      )      
    }
  }, ignoreInit = TRUE)
    
  create_plot <- eventReactive(list(input$refreshPlot, input$log), {
    series <- myseries()
    req(series)
    
    if (series %>% filter(!is.na(series)) %>% nrow()) {
      mydata <- parse_data() %>%
        filter(channel == input$channel)
      
      c_lim <- try(cellranger::as.cell_limits(input$range), silent=TRUE)
      
      if (!inherits(c_lim,"try-error"))
        mydata <- semi_join(mydata, range_to_wells(c_lim),
                            by=c("row","col"))
      
      if (input$blank > 0)
        mydata <- mutate(mydata, value=value-input$blank)
      
      pl <- mydata %>%
        group_by(series,row,col) %>%
        mutate(time_sec=as.numeric(datetime - min(datetime))) %>%
        # identity() %>% print()
        ggplot() +
        facet_grid(row~col) +
        geom_point(aes(time_sec/3600,value,col=series),
                   size=1, stroke=0, alpha=.7) +
        labs(x="time (h)") +
        theme_bw() +
        theme(legend.position='bottom') +
        NULL
      
      if (input$log)
        pl <- pl + scale_y_log10()
      
      return(pl)
    }
  })
  
  output$gc_plot <- plotly::renderPlotly({    
    pl <- create_plot()
    req(pl)
    plotly::ggplotly(pl, height=720) %>%
      plotly::layout(legend=list(orientation="h", x=0, y=-0.07, itemsizing="constant"))
  })
}

# RUN APP ####
shinyApp(ui, server)
