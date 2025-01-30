library(shiny)

source('ui.R')
source('server.R')

# initialize logger
source("logger.R")
init_logger()

shinyApp(ui = ui, server = server)