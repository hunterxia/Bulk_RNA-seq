#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#
library(shiny)

server <- function(input, output, session) {
  dataset <- mainTabServer("main")
  qcTabServer("qc", dataset)
  PCACorrelationTabServer("pca_correlation", dataset)
  clusteringTabServer("clustering")
  PairwiseComparisonTabServer("pariwise_comparison",dataset)
}