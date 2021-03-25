# This is a Shiny web application. You can run the application by executing `vngGrowthCurves:::Momentum_GC_run()`

library(shiny)
# library(shinyFiles)
library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)
library(vngGrowthCurves)

# Define server logic required to parse files and plot the growth curves
# use a function factory to customise the paths in the UI
shinyServer(
  function(input, output, session) {

  # IMPORTANT! close session when running standalone
  # this is needed to terminate the R process when the
  # shiny app session ends. Otherwise, you end up with a zombie process
  if (!interactive()) {
    session$onSessionEnded(function() {
      stopApp()
      q("no")
    })
  }
  
  # enable directory selection in UI
  volumes <- c(Home = fs::path_home(), 
               # "R Installation" = R.home(), 
               shinyFiles::getVolumes()() )
  dir <- getOption("MomentumGCdir")
  if (!is.null(dir)) volumes <- c(rlang::set_names(fs::path_tidy(dir), fs::path_file(dir)), 
                                  volumes)
  shinyFiles::shinyDirChoose(
    input,
    'path',
    roots = volumes,
    session=session
  )
  output$path <- renderPrint({
    if (is.integer(input$path)) {
      cat("No input directory has been selected.")
    } else {
      shinyFiles::parseDirPath(volumes, input$path)
    }
  })
  
  myseries <- reactive({
  # myseries is NULL when no file can be parsed
    if (is.integer(input$path)) {
      NULL
    } else {
      re <- if (input$re=="") NULL else input$re
      # browser()
      .paths <- shinyFiles::parseDirPath(volumes, input$path) %>% 
        { try(fs::dir_ls(., type="file", regexp=re)) } #dir_ls() crashes when regex is not correct
      
      if (class(.paths) == "try-error") NULL else
        tibble(path=.paths) %>%
        mutate(filename = fs::path_file(path),
               is_txt=stringr::str_detect(filename, "\\.txt$"),
               series = stringr::str_extract(filename, "^.*_pl\\d+")
        ) %>%
        group_by(series) %>%
        nest() %>% 
        mutate(channels=map(data, ~try( read_Biotek_Synergy2_matrices(.$path[.$is_txt][1], .ch_only=T) ))) %>%
        identity()
    }
  })
  
  output$series <- renderPrint({
    if (is.null(myseries())) {
      invisible() # return nothing (as per the doc)
    } else if (myseries() %>% pull(channels) %>% unique() %>% length() > 1) {
      # case of channels mismatch
      print("Error: found files with different measurement channels, please adjust your selection criteria.")
    } else {
      # case of parsable data
      if (any(is.na(myseries()$series)))
        print("Message: found files not matching the pattern...")
      
      if (myseries() %>% filter(is.na(series)) %>% nrow())
        if (myseries() %>% filter(is.na(series)) %>% unnest(data) %>% filter(!is_txt) %>% nrow()) 
          print("Warning: ignoring files with extension different from `.txt`...")

      cat("Found ", myseries() %>% filter(!is.na(series)) %>% nrow(), " files series:")
      if (myseries() %>% filter(!is.na(series)) %>% nrow())
        myseries() %>% filter(!is.na(series)) %>% 
        unnest(data) %>% filter(is_txt) %>% 
        group_by(series) %>% summarise(n_files = n(), .groups="drop_last") %>%
        knitr::kable(format='html')
    }
  })
  
  range_to_wells <- function(c_lim) {
    # NB: excel cells use letters for columns while microplate wells use them for rows... 
    expand_grid(col=c_lim$ul[1]:c_lim$lr[1], row=LETTERS[c_lim$ul[2]:c_lim$lr[2]])
  }
  
  parse_data <- reactive({
    if (!is.null(myseries())) # redundant with is_empty()
      if (!is_empty(myseries()))
        if (myseries() %>% filter(!is.na(series)) %>% nrow()) {
          myseries() %>% filter(!is.na(series)) %>% select(-channels) %>% 
            unnest(data) %>% filter(is_txt) %>% 
            extract(filename, c('filedate', 'filetime'), "_pl\\d+_(\\d{8})_(\\d{6})(?:_FE_)?.txt$", remove=FALSE) %>% 
            rowwise() %>% 
            mutate(
              datetime=lubridate::ymd_hms(paste0(filedate, filetime)), filedate=NULL, filetime=NULL,
              data=list(vngGrowthCurves::read_Biotek_Synergy2_matrices(path)),
              # data=purrr::map(path, ~vngGrowthCurves::read_Biotek_Synergy2_matrix(.))
            ) %>%
            unnest(data) %>%
            unnest(data) %>%
            mutate(value=as.numeric(value)) %>% 
            identity()
        }
  })
  
  observe({
    if (!is.null(myseries())) {
      .ch <- myseries() %>% pull(channels) %>% unique()
      
      if (length(.ch) == 1 && typeof(.ch[[1]])=="character")
        updateSelectInput(session, "channel",
                          label = "Channel",
                          choices = .ch[[1]],
                          selected = .ch[[1]][1]
        )
    }
  })
        
  create_plot <- eventReactive(c(input$refreshPlot, input$log), {
    if (!is.null(myseries()))
      if (myseries() %>% filter(!is.na(series)) %>% nrow()) {
        mydata <- parse_data() %>% 
          filter(channel == input$channel)
        
        c_lim <- try(cellranger::as.cell_limits(input$range), silent = TRUE)
        if (class(c_lim)[1] != 'try-error')
          mydata <- semi_join(mydata, range_to_wells(c_lim), by=c("row", "col"))
        
        if (input$blank > 0)
          mydata <- mutate(mydata, value=value-input$blank)
          
        pl <- mydata %>% 
          group_by(series, row, col) %>%
          mutate(time_sec=as.numeric(datetime - min(datetime))) %>%
          # identity() %>% print()
          ggplot() +
          facet_grid(row~col) +
          geom_point(aes(time_sec / 3600, value, col=series), size=1, stroke=0, alpha=.7) +
          labs(x="time (h)") +
          theme_bw() +
          theme(legend.position = 'bottom') +
          NULL
        
        if (input$log)
          pl <- pl + scale_y_log10()
        
        return(pl)
      }
  })
  
  output$gc_plot <- plotly::renderPlotly({
    if (!is.null(create_plot()))
        plotly::ggplotly(create_plot(), height=720) %>%
          plotly::layout(legend = list(orientation = "h", x=0, y=-0.07, itemsizing="constant"))
  })  
  
})
