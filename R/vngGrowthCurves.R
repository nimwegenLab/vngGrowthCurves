.onLoad <- function(libname, pkgname) {
# put here was you need to execute at package loading
  
  invisible()
}

globalVariables(c(".", ":=", "time", "hours", "sec", "channel", "value", "well", "last_col", "step", "id_first", "id_last",
                  "a", "i", "sd", "id_line", "Time"))

create_empty_plate <- function(.nrow=8, .ncol=12, ...) {
  dplyr::left_join(
    dplyr::tibble(a=1, ..., row=LETTERS[1:.nrow]), 
    dplyr::tibble(a=1, col=1:.ncol), 
    by="a") %>%
    dplyr::select(-a) %>% 
    dplyr::mutate(well=paste0(row, col))
}

add_col_var <- function(.df, .name, .values, .row_subset=NULL) {
  .ncol <- length(unique(.df$col))
  
# .row_subset is a vector of row indices, e.g. `LETTERS[1:4]`
  if (length(.values) < .ncol) warning("`.values` has less than ", .ncol, " values (values aren't recycled).")
  if (length(.values) > .ncol) warning("`.values` has more than ", .ncol, " values (values beyond the first twelve are discarded).")
  
  if (!is.null(.row_subset)) {
    if (length(setdiff(.row_subset, LETTERS[1:8])) > 0) warning("`.row_subset` has values not matching row indices (they will be ignored).")
    .row_subset <- intersect(.row_subset, LETTERS[1:8])
  } else {
    # better surcharge the cross_join table otherwise the columns by which to left_join are not always the same
    .row_subset <- LETTERS[1:8]
  }
  
  dplyr::bind_rows(
    dplyr::filter(.df, ! row %in% .row_subset),
    dplyr::left_join(
      dplyr::filter(.df, row %in% .row_subset) %>% 
        dplyr::mutate({{.name}} := NULL), # override existing values
      dplyr::tibble(col=1:length(.values)) %>% 
        dplyr::mutate({{.name}} := .values) %>% 
        dplyr::cross_join(dplyr::tibble(row=.row_subset)),
      by = c("row", "col"),
    )
  ) 
  # TODO: sort by id (e.g. date), row, col (figure out how to identify which column was used as id)
}

add_row_var <- function(.df, .name, .values, .col_subset=NULL) {
  .nrow <- length(unique(.df$row))
  
# .col_subset is a vector of column indices, e.g. `1:4`
  if (length(.values) < .nrow) warning("`.values` has less than ", .nrow, " values (values aren't recycled).")
  if (length(.values) > .nrow) warning("`.values` has more than ", .nrow, " values (values beyond the first eight are discarded).")
  
  if (!is.null(.col_subset)) {
    if (length(setdiff(.col_subset, 1:12)) > 0) warning("`.col_subset` has values not matching column indices (they will be ignored).")
    .col_subset <- intersect(.col_subset, 1:12)
  } else {
    # better surcharge the cross_join table otherwise the columns by whihc to left_join are not always the same
    .col_subset <- 1:12
  }

  # browser()
  # if (tibble::has_name(.df, {{.name}}) ) # did not find the right way to handle `.name` for this
  #   if (dplyr::filter(.df, col %in% .col_subset) %>% dplyr::filter(!is.na({{.name}})) %>% nrow() > 0)
  #     warning("Values are overriden by add_row_var()." )
  
  dplyr::bind_rows(
    dplyr::filter(.df, ! col %in% .col_subset),
    dplyr::left_join(
      dplyr::filter(.df, col %in% .col_subset) %>% 
        dplyr::mutate({{.name}} := NULL), # override existing values
      dplyr::tibble(row=LETTERS[1:length(.values)]) %>% 
        dplyr::mutate({{.name}} := .values) %>% 
        dplyr::cross_join(dplyr::tibble(col=.col_subset)),
      by = c("row", "col"),
    )
  ) 
  # TODO: sort by id (e.g. date), row, col (figure out how to identify which column was used as id)
}

read_Biotek_Synergy2_matrices <- function(.path, .channels = "all", .ch_only=FALSE) {
  # .channels is either a vector of channels to be kept (one-based indices, or names), or "all"
  
  warning("read_Biotek_Synergy2_matrices() is deprecated; use read_Biotek_Synergy2_kinetic() instead.")
  .lines <- readLines(.path)

  .ids <- stringr::str_which(.lines, "^$") # find empty lines
  if (.ch_only) 
    return(.lines[.ids - 10])
  if (.channels != "all") .ids <- .ids[.channels]
  
  # identify delimiter (tab and comma supported)
  .delim <- ""
  if (stringr::str_count(.lines[4], "\\t") == 13) { 
    .delim <- "\t" 
    } else if (stringr::str_count(.lines[4], ",") == 13) { .delim <- "," }
  if (.delim == "") stop("text delimiter not recognized")
  
  parse_table_text <- function(text, delim)
    utils::read.table(text=text, sep=delim, header=FALSE, stringsAsFactors=FALSE) %>% 
    # stats::setNames(c("row", 1:12, 'last_col')) %>% 
    # dplyr::select(-last_col) %>% 
    # handle cases where the last column might have separator chars
    dplyr::select(1:13) %>%
    stats::setNames(c("row", 1:12)) %>%
    # reshape wide to long
    tidyr::gather(col, value, dplyr::matches("\\d+")) %>% 
    dplyr::mutate(well=paste0(row, col), col=as.numeric(col))
  
  dplyr::tibble(id_line = .ids - 8) %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(channel = purrr::map_chr(id_line, ~.lines[.-2])) %>% 
    dplyr::mutate(data = purrr::map(id_line, ~ {
      paste(.lines[seq(.x, .x+7)], collapse='\n') %>%
        parse_table_text(delim=.delim)
    }),
    id_line = NULL,
    )
}

read_Biotek_Synergy2_matrix <- function(.path) {
  # read only the first channel, ignore silently following ones
  read_Biotek_Synergy2_matrices(.path, .channels=1) %>% 
    tidyr::unnest(data)
}

read_Biotek_Synergy2_kinetic <- function(.path, ...) {
#  read_spec_kinetic() is written such as to call read.table only once 
# (much faster than calling it for each timepoint and channel)
  .lines <- readLines(.path)
  .l_idx <- stringr::str_detect(.lines, "(?:Kinetic read|Time) (\\d+) (.*)") %>% which %>% (function(.x) .x-1) 

  .data <- lapply(.l_idx, function(.i) {
    # extract the channel name
    .m <- stringr::str_match(.lines[.i], "(?:Read \\d+:)*(.+)")
    .ch <- .m[2]
    # extract timelapse step and timestamp
    .m <- stringr::str_match(.lines[.i+1], "(?:Kinetic read|Time) (\\d+) \\((.*)\\)")
    .time <- .m[3]
    .step <- .m[2] # add time step (because different channels are recorded at different times)
    # prepend to rows and concatenate
    paste(.ch, .time, .step, .lines[(.i+3):(.i+10)], sep='\t') %>% 
      paste(collapse='\n')
  }) %>% 
    paste(collapse='\n') %>% 
    stringr::str_replace_all(stringr::fixed("OVRFLW"), "Inf") %>% 
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
    dplyr::mutate(well=paste0(row, col), col=as.integer(col), channel=as.character(channel))
  
 if(nrow(dplyr::filter(.data, step>1, time==0))) 
   warning('CRITICAL: You should check that your data file doesnt contain empty measurements at its end (this happens when the acquisition is stopped manually...). It is strongly advised to delete those manually and import again.')
 
  return(
    .data %>% 
      dplyr::filter(time > 0) %>% 
      dplyr::full_join(dplyr::tibble(...), by=character())
    )
}

read_Biotek_Synergy2_columns <- function(.path, .encoding="latin1", .skip=0, ...) {
  # browser()
  .lines <- readr::read_lines(.path, locale=readr::locale(encoding=.encoding), skip=.skip)
  .ids <- stringr::str_which(.lines, "^[^\\t]+$") # find empty lines (not containing `TAB`)
  
  parse_table_text <- function(.text)
    stringr::str_replace_all(.text, stringr::fixed("OVRFLW"), "Inf") %>% 
    readr::read_delim(
      delim='\t', 
      col_types=readr::cols(readr::col_time("%h:%M:%S"), .default=readr::col_double())
    ) %>% 
    dplyr::mutate(time=as.numeric(Time), step=dplyr::row_number()) %>% 
    dplyr::select(-Time, -dplyr::contains("T\u00B0")) %>% 
    dplyr::filter(!(dplyr::row_number() > 1 & time==0)) %>% 
    tidyr::pivot_longer(cols=c(-time, -step), names_to="well", values_to = "value")

  dplyr::tibble(id_first = .ids + 2, 
                id_last = dplyr::lead(.ids+-1, default=length(.lines)) ) %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(channel = purrr::map_chr(id_first, ~.lines[.-2])) %>% 
    dplyr::mutate(
      data = purrr::map2(id_first, id_last, ~ {
        paste(.lines[.x:.y], collapse='\n') %>% 
          parse_table_text()
      }),
      id_first = NULL, id_last=NULL,
    ) %>% 
    tidyr::unnest(data) %>% 
    tidyr::extract(well, c("row", "col"), "([A-H])(\\d+)", remove=FALSE, convert=TRUE) %>% 
    dplyr::full_join(dplyr::tibble(...), by=character())
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

time_at_od <- function(.time, .od, .target_od) {
  .idx <- which(.od > .target_od & lag(.od) < .target_od)
  if (length(.idx) > 0) return(.time[max(.idx)]) else return(max(.time))
}


MomentumGC_run <- function(dir=NULL, ...)
  shiny::runApp(system.file("MomentumGC", package="vngGrowthCurves"), ...)
