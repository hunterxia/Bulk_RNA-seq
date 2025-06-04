#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

# ui.R

library(shiny)
source("./R/main.R")
source("./R/qc.R")
source("./R/pca_correlation.R")
source("./R/pairwise_comparison.R")
source("./R/clustering.R")
source("./R/individual_gene.R")
ui <- fluidPage(
  titlePanel("Bulk RNA-seq Analysis App"),
  tabsetPanel(
    tabPanel(
      "Main",
      mainTabUI("main")
    ),
    tabPanel(
     "Quality Control",
     qcTabUI("qc")
    ),
    tabPanel(
      "Global Comparison",
      PCACorrelationTabUI("pca_correlation")
    ),
    tabPanel(
      "Pairwise Comparison",
      PairwiseComparisonTabUI("pairwise_comparison")
    ),
    tabPanel(
      "Clustering",
      clusteringTabUI("clustering")
    ),
    tabPanel("Individual Gene", individualGeneTabUI("individual_gene"))
  )
)
