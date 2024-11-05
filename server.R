#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#
library(shiny)
source("./R/main.R")
source("./R/qc.R")
source("./R/pca_correlation.R")
source("./R/pairwise _comparison.R")
source("./R/clustering.R")
server <- function(input, output, session) {
  session$onSessionEnded(function() {
    while (!is.null(dev.list())) {
      dev.off()
      graphics.off() 
    }
  })
  dataset <- mainTabServer("main")
  qcTabServer("qc", dataset)
  PCACorrelationTabServer("pca_correlation", dataset)
  clustering_dataset <- clusteringTabServer("clustering", dataset)
  PairwiseComparisonTabServer("pairwise_comparison", dataset)
}