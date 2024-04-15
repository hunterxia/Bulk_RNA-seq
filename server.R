#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#
library(shiny)
source("R/main.R")

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  dataset <- mainTabServer("mainTabModule")
  qcTabServer("qc", dataset)
  PCACorrelationTabServer("pca_correlation")
  clusteringTabServer("clustering")
  PairwiseComparisonTabServer("pariwise_comparison")
}