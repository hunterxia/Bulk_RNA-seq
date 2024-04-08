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
library(DT) 

# Define UI for main tab
mainTabUI <- function() {
  tabsetPanel(
    tabPanel("Main", 
             fluidPage(
               titlePanel("Bulk RNA-seq"),
               
               dataTableOutput("expression_table"),
               
               dataTableOutput("groups_table"),
               
               shiny::checkboxInput("toggle_groups", "Toggle Groups On/Off", value = FALSE),
              
               uiOutput("sample_checkboxes"),
              
               fluidRow(
                 column(3, numericInput("cutoff_expression", "Cutoff Expression", value = 0)),
                 column(3, numericInput("samples_cutoff", "Samples Cutoff", value = 0))
               ),
               selectInput("dataset", label = "Dataset", choices = ls("package:datasets")),
               verbatimTextOutput("summary"),
               tableOutput("table")
             )
    ),
    
    tabPanel("QC Tab",
             fluidPage(
               
             )
    ),
    
    tabPanel("PCA and Correlation",
             fluidPage(
               
             )
    ),
    tabPanel("Pairwise Comparison",
             PairwiseComparisonTabUI("pariwise_comparison")
    ),
    tabPanel("Clustering",
             clusteringTabUI("clustering")
    )
   
  )
}
