.onLoad <- function(libname, pkgname) {
# put here was you need to execute at package loading
  
  invisible()
}

globalVariables(c(".", "time", "hours", "sec", "channel", "value", "last_col", "step",
                  "i", "sd"))

read_Biotek_Synergy2_matrix <- function(.path) {
  .lines <- readLines(.path)
  # extract the channel name
  paste(.lines[3:length(.lines)], collapse='\n') %>% 
    utils::read.table(text=., sep="\t", header=FALSE, stringsAsFactors=FALSE) %>% 
    stats::setNames(c("row", 1:12, 'last_col')) %>% 
    dplyr::select(-last_col) %>% 
    # reshape wide to long
    tidyr::gather(col, value, dplyr::matches("\\d+")) %>% 
    dplyr::mutate(well=paste0(row, col), col=as.numeric(col), channel=.lines[1])
}

read_Biotek_Synergy2_kinetic <- function(.path) {
#  read_spec_kinetic() is written such as to call read.table only once 
# (much faster than calling it for each timepoint and channel)
  .lines <- readLines(.path)
  .l_idx <- stringr::str_detect(.lines, "(?:Kinetic read)|(?:Time) (\\d+) (.*)") %>% which %>% (function(.x) .x-1)

  .data <- lapply(.l_idx, function(.i) {
    # extract the channel name
    .m <- stringr::str_match(.lines[.i], "(?:Read \\d+:)*(.+)")
    .ch <- .m[2]
    # extract timelapse step and timestamp
    .m <- stringr::str_match(.lines[.i+1], "(?:Kinetic read)|(?:Time) (\\d+) \\((.*)\\)")
    .time <- .m[3]
    .step <- .m[2] # add time step (because different channels are recorded at different times)
    # prepend to rows and concatenate
    paste(.ch, .time, .step, .lines[(.i+3):(.i+10)], sep='\t') %>% 
      paste(collapse='\n')
  }) %>% 
    paste(collapse='\n') %>% 
    utils::read.table(text=., sep="\t", header=FALSE, stringsAsFactors=FALSE) %>% 
    stats::setNames(c("channel", "time", "step", "row", 1:12, 'last_col')) %>% 
    dplyr::select(-last_col) %>% 
    # convert time to float
    # dplyr::mutate(time=lubridate::hms(time)) %>%
    tidyr::extract(time, c('hours', 'min', 'sec'), '(\\d+):(\\d+):(\\d+)') %>%
    dplyr::mutate(time=as.numeric(hours)*3600 + as.numeric(min)*60 + as.numeric(sec)) %>%
    dplyr::select(-hours, -min, -sec) %>%
    # reshape wide to long
    tidyr::gather(col, value, dplyr::matches("\\d+")) %>% 
    dplyr::mutate(well=paste0(row, col), col=as.numeric(col), channel=as.character(channel))
  
 if(nrow(dplyr::filter(.data, step>1, time==0))) 
   warning('CRITICAL: You should check that your data file doesnt contain empty measurements at its end (this happens when the acquisition is stopped manually...). It is strongly advised to delete those manually an import again.')
 
  return(.data)
}


find_blank_od <- function(.t, .od, .tmin=0, .tmax=Inf, .npoints=10, .nstep=5, .nthresh=50, .cv_thresh=.01) {
  stopifnot(is.numeric(.t),
            is.numeric(.od),
            length(.t) == length(.od))
  # browser()
  .stats <- data.frame(i=seq(which.min(abs(.t-.tmin)), which.min(abs(.t-.tmax)), .nstep)) %>%
    dplyr::mutate(t=.t[i]) %>% 
    dplyr::group_by(i, t) %>% 
    dplyr::do((function(.df){
      .tmp <- .od[.df$i:(.df$i+.npoints)]
      data.frame(mean=mean(.tmp), sd=stats::sd(.tmp))
      })(.))
  
  .out <- .stats %>% dplyr::ungroup() %>%
    dplyr::filter(sd/mean < .cv_thresh) %>% 
    dplyr::filter(mean==min(mean)) %>% 
    dplyr::slice(1) %>% # if several are found use the first one
    with(list(value=mean, time=t))
  
  .n_below <- sum(.od<.out$blank)
  if (.n_below > .nthresh || length(.out$value) != 1) return(list(value=NA, time=NA))
  
  # todo: check that od is stable for a certain time after time_blank
  .out
}
