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

ui <- fluidPage(
  titlePanel("RNA-seq Data Visualization App"),
  tabsetPanel(
    tabPanel("MainTab",
             mainTabUI("main")),
    tabPanel("QC Tab",
             qcTabUI("qc")),
    tabPanel("PCA and Correlation",
             PCACorrelationTabUI("pca_correlation")),
    tabPanel("Pairwise Comparison",
             PairwiseComparisonTabUI("pairwise_comparison")),
    tabPanel("Clustering",
             clusteringTabUI("clustering"))
  )
)
