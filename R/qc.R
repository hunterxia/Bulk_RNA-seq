library(viridis)
library(hrbrthemes)
library(ggplot2)
library(tidyr)

qcTabServer <- function(id, dataset) {
  moduleServer(id, function(input, output, session) {
    create_hist_plot <- function() {
      df <- dataset$expression_data()
      df_log_normalized <- log1p(df[3:ncol(df)])
      df_melt <- melt(df_log_normalized, variable.factor=TRUE, variable.name="variable", value.name="value")
      bins_pt5 <- table(cut(df_melt$value, breaks=seq(0, 16.5, by=input$bin_size)), df_melt$variable)
      bins_pt5 <- as.data.frame(bins_pt5)
      bins_pt5 %>%
        ggplot(aes(x=Var1, y=Freq, group=Var2, color=Var2)) +
        geom_line() +
        scale_color_viridis(discrete = TRUE) +
        ggtitle("Popularity of American names in the previous 30 years") +
        theme_ipsum() + 
        labs(
          x = paste0("Bins(binsize=", input$bin_size, ")"),
          y = "Gene Counts",
          title = "Histogram before filtering"
        ) +
        guides(color = guide_legend(nrow = 6, byrow = TRUE, title.position = "top")) +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title.x = element_text(vjust = 0.5, hjust = 0.5),
          axis.title.y = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          legend.box.background = element_rect(colour = "black", fill = NA, size = 0.2),
          legend.title = element_blank(),
          legend.position = c(1, 1),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
          legend.key.size = unit(0.5, "cm")
        )
    }
    
    create_filtered_hist_plot <- function() {
      df <- dataset$filtered_data()
      df_log_normalized <- log1p(df[3:ncol(df)])
      df_melt <- melt(df_log_normalized, variable.factor=TRUE, variable.name="variable", value.name="value")
      bins_pt5 <- table(cut(df_melt$value, breaks=seq(0, 16.5, by=input$bin_size)), df_melt$variable)
      bins_pt5 <- as.data.frame(bins_pt5)
      bins_pt5 %>%
        ggplot(aes(x=Var1, y=Freq, group=Var2, color=Var2)) +
        geom_line() +
        scale_color_viridis(discrete = TRUE) +
        ggtitle("Popularity of American names in the previous 30 years") +
        theme_ipsum() + 
        labs(
          x = paste0("Bins(binsize=", input$bin_size, ")"),
          y = "Gene Counts",
          title = "Histogram after filtering"
        ) +
        guides(color = guide_legend(nrow = 6, byrow = TRUE, title.position = "top")) +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title.x = element_text(vjust = 0.5, hjust = 0.5),
          axis.title.y = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          legend.box.background = element_rect(colour = "black", fill = NA, size = 0.2),
          legend.title = element_blank(),
          legend.position = c(1, 1),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
          legend.key.size = unit(0.5, "cm")
        )
    }
    
    create_grouped_plot <- function() {
      df <- dataset$expression_data()
      groups <- dataset$groups_data()
      df_long <- pivot_longer(df, cols = -names(df)[1:2], names_to = "Variables", values_to = "Value")
      df_long <- df_long %>%
        left_join(groups, by = setNames( names(groups)[1], "Variables")) 
      df_long <- df_long %>%
        group_by(df_long[,1], df_long[,2], df_long[,5]) %>%
        summarize(Value = sum(Value), .groups = "drop")
      result_data <- pivot_wider(df_long, names_from = names(df_long)[3], values_from = Value, values_fill = list(value = 0))
      
      df_melt <- melt(result_data, variable.factor=TRUE, variable.name="variable", value.name="value")
      bins_pt5 <- table(cut(df_melt$value, breaks=seq(0, 16.5, by=input$bin_size)), df_melt$variable)
      bins_pt5 <- as.data.frame(bins_pt5)
      
      bins_pt5 %>%
        ggplot(aes(x=Var1, y=Freq, group=Var2, color=Var2)) +
        geom_line() +
        scale_color_viridis(discrete = TRUE) +
        ggtitle("Popularity of American names in the previous 30 years") +
        theme_ipsum() + 
        labs(
          x = paste0("Bins(binsize=", input$bin_size, ")"),
          y = "Gene Counts",
          title = "Histogram after grouping"
        ) +
        guides(color = guide_legend(nrow = 6, byrow = TRUE, title.position = "top")) +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title.x = element_text(vjust = 0.5, hjust = 0.5),
          axis.title.y = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          legend.box.background = element_rect(colour = "black", fill = NA, size = 0.2),
          legend.title = element_blank(),
          legend.position = c(1, 1),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
          legend.key.size = unit(0.5, "cm") 
        )
    }
    
    create_grouped_filtered_plot <- function() {
      df <- dataset$filtered_data()
      groups <- dataset$groups_data()
      df_long <- pivot_longer(df, cols = -names(df)[1:2], names_to = "Variables", values_to = "Value")
      df_long <- df_long %>%
        left_join(groups, by = setNames( names(groups)[1], "Variables")) 
      df_long <- df_long %>%
        group_by(df_long[,1], df_long[,2], df_long[,5]) %>%
        summarize(Value = sum(Value), .groups = "drop")
      
      result_data <- pivot_wider(df_long, names_from = names(df_long)[3], values_from = Value, values_fill = list(value = 0))
      
      df_melt <- melt(result_data, variable.factor=TRUE, variable.name="variable", value.name="value")
      bins_pt5 <- table(cut(df_melt$value, breaks=seq(0, 16.5, by=input$bin_size)), df_melt$variable)
      bins_pt5 <- as.data.frame(bins_pt5)
      
      bins_pt5 %>%
        ggplot(aes(x=Var1, y=Freq, group=Var2, color=Var2)) +
        geom_line() +
        scale_color_viridis(discrete = TRUE) +
        ggtitle("Popularity of American names in the previous 30 years") +
        theme_ipsum() + 
        labs(
          x = paste0("Bins(binsize=", input$bin_size, ")"),
          y = "Gene Counts", 
          title = "Histogram after grouping and filtering"
        ) +
        guides(color = guide_legend(nrow = 6, byrow = TRUE, title.position = "top")) +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title.x = element_text(vjust = 0.5, hjust = 0.5), 
          axis.title.y = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          legend.box.background = element_rect(colour = "black", fill = NA, size = 0.2),
          legend.title = element_blank(),
          legend.position = c(1, 1),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
          legend.key.size = unit(0.5, "cm")
        )
      
    }
    
    output$hist <- renderPlot({
      req(dataset$expression_data())
      create_hist_plot()
    }, res = 96)
    
    output$filtered_hist <- renderPlot({
      req(dataset$expression_data())
      create_filtered_hist_plot()
    }, res = 96)
    
    output$grouped_hist <- renderPlot({
      req(dataset$expression_data())
      create_grouped_plot()
    }, res = 96)
    
    
    output$grouped_filtered_hist <- renderPlot({
      req(dataset$expression_data())
      create_grouped_filtered_plot()
    }, res = 96)
  
    
    output$download_hist <- downloadHandler(
      filename = function() {
        paste("plots-", Sys.Date(), ".zip", sep = "")
      },
      content = function(file) {
        tempDir <- "www"
        dir.create(tempDir, recursive = TRUE, showWarnings = FALSE)
        filePath1 <- file.path(tempDir, "plot1.png")
        filePath2 <- file.path(tempDir, "plot2.png")
        filePath3 <- file.path(tempDir, "plot3.png")
        filePath4 <- file.path(tempDir, "plot4.png")
        # Save plots
        ggsave(filePath1, plot = create_hist_plot(), device = "png", width = 10, height = 8, dpi = 300)
        ggsave(filePath2, plot = create_filtered_hist_plot(), device = "png", width = 10, height = 8, dpi = 300)
        ggsave(filePath3, plot = create_grouped_plot(), device = "png", width = 10, height = 8, dpi = 300)
        ggsave(filePath4, plot = create_grouped_filtered_plot(), device = "png", width = 10, height = 8, dpi = 300)
        # Create zip file
        zip(file, files = c(filePath1, filePath2, filePath3, filePath4))
      }
    )
    
  })
}

qcTabUI <- function(id) {
  fluidPage(
    titlePanel("QC"),
    fluidRow(
      column(3, 
        sliderInput(NS(id, "bin_size"), "Bin size", min = 0, max = 1, value = 0.25, step = 0.25),
        downloadButton(NS(id, "download_hist"), "Download Plots")
      ),
      column(3, 
        
      ),
    ),
    fluidRow(
      column(6, 
        plotOutput(NS(id, "hist"))
      ),
      column(6, 
        plotOutput(NS(id, "filtered_hist"))
      )
    ),
    fluidRow(
      column(6, 
        plotOutput(NS(id, "grouped_hist"))
      ),
      column(6, 
        plotOutput(NS(id, "grouped_filtered_hist"))
      )
    )
  )
}