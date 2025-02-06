library(shiny)
library(logger)

source('ui.R')
source('server.R')

shinyApp(ui = ui, server = server)