library(logger)

init_logger <- function(log_dir = "logs") {
  # set log directory
  if (missing(log_dir)) {
    log_dir <- file.path(getwd(), "logs")
  }

  if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # set up log file
  timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
  log_file <- file.path(log_dir, paste0("shiny_", timestamp, ".log"))

  if (!file.exists(log_dir)) {
    warning("Cannot create log directory, using console output")
    log_appender(appender_console)
  } else {
    # set log output to file
    log_appender(appender_file(log_file))
  }

  # record initial log
  try(log_info("Shiny session started"))
}