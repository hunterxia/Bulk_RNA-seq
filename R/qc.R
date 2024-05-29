library(viridis)
library(hrbrthemes)
library(ggplot2)
library(tidyr)
library(shinyWidgets)
library(plotly)

qcTabServer <- function(id, dataset) {
  moduleServer(id, function(input, output, session) {
    create_hist_plot <- function() {
      df <- dataset$expression_data()
      groups <- dataset$groups_data()
      selected_samples <- dataset$selected_samples_data()
      
      df <- df %>%
        select(c('Symbols', 'Genes', all_of(selected_samples)))
      
      colnames(df)[1] <- "Symbols"
      colnames(df)[2] <- "Genes"
      colnames(groups)[1] <- "Samples"
      colnames(groups)[2] <- "EXPERIMENTAL_GROUP"
      
      
      if (input$grouped) {
        # show by grouped
        df_long <- pivot_longer(df, cols = -names(df)[1:2], names_to = "Samples", values_to = "Values")
        
        df_long_grouped <- df_long %>%
          left_join(groups, by = "Samples") 
        
        df_long_grouped_sum <- df_long_grouped %>%
          group_by(Genes, Symbols, EXPERIMENTAL_GROUP) %>%
          summarize(Values = sum(Values), .groups = "drop")
        
        df_grouped <- pivot_wider(df_long_grouped_sum, names_from = EXPERIMENTAL_GROUP, values_from = Values, values_fill = list(value = 0))
        
        df_log_normalized <- log1p(df_grouped[3:ncol(df_grouped)])
        df_melt <- melt(df_log_normalized, variable.factor=TRUE, variable.name="variable", value.name="value")
        bins_pt5 <- table(cut(df_melt$value, breaks=seq(0, 16.5, by=input$bin_size)), df_melt$variable)
        bins_pt5 <- as.data.frame(bins_pt5)
        bins_pt5 <- bins_pt5 %>%
          left_join(groups, by = c("Var2" = "EXPERIMENTAL_GROUP")) 
        
      } else {
        # show by smaple
        df_log_normalized <- log1p(df[3:ncol(df)])
        df_melt <- melt(df_log_normalized, variable.factor=TRUE, variable.name="variable", value.name="value")
        bins_pt5 <- table(cut(df_melt$value, breaks=seq(0, 16.5, by=input$bin_size)), df_melt$variable)
        bins_pt5 <- as.data.frame(bins_pt5)
      }

      if("Color" %in% colnames(bins_pt5)) {
        colors <- unique(bins_pt5$Color)
        names(colors) <- unique(bins_pt5$Var2)
        scale_color <- scale_color_manual(values = colors)
      } else {
        scale_color <- scale_color_viridis(discrete = TRUE)
      }
      
      p <- bins_pt5 %>%
        ggplot(aes(x=Var1, y=Freq, group=Var2, color=Var2, text = paste(ifelse(input$grouped, "Group:", "Sample:"), Var2, "<br>Interval:", Var1, "<br>Count:", Freq))) +
        geom_line() +
        scale_color +
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
      
        ggplotly(p, tooltip = "text")
    }
    
    create_filtered_hist_plot <- function() {
      df <- dataset$filtered_data()
      groups <- dataset$groups_data()
      
      colnames(df)[1] <- "Symbols"
      colnames(df)[2] <- "Genes"
      colnames(groups)[1] <- "Samples"
      colnames(groups)[2] <- "EXPERIMENTAL_GROUP"
      
      if (input$grouped) {
        # show by grouped
        df_long <- pivot_longer(df, cols = -names(df)[1:2], names_to = "Samples", values_to = "Values")
        
        df_long_grouped <- df_long %>%
          left_join(groups, by = "Samples") 
        
        df_long_grouped_sum <- df_long_grouped %>%
          group_by(Genes, Symbols, EXPERIMENTAL_GROUP) %>%
          summarize(Values = sum(Values), .groups = "drop")
        
        df_grouped <- pivot_wider(df_long_grouped_sum, names_from = EXPERIMENTAL_GROUP, values_from = Values, values_fill = list(value = 0))
        
        df_log_normalized <- log1p(df_grouped[3:ncol(df_grouped)])
        df_melt <- melt(df_log_normalized, variable.factor=TRUE, variable.name="variable", value.name="value")
        bins_pt5 <- table(cut(df_melt$value, breaks=seq(0, 16.5, by=input$bin_size)), df_melt$variable)
        bins_pt5 <- as.data.frame(bins_pt5)
        bins_pt5 <- bins_pt5 %>%
          left_join(groups, by = c("Var2" = "EXPERIMENTAL_GROUP")) 
        
      } else {
        # show by smaple
        df_log_normalized <- log1p(df[3:ncol(df)])
        df_melt <- melt(df_log_normalized, variable.factor=TRUE, variable.name="variable", value.name="value")
        bins_pt5 <- table(cut(df_melt$value, breaks=seq(0, 16.5, by=input$bin_size)), df_melt$variable)
        bins_pt5 <- as.data.frame(bins_pt5)
      }
      
      if("Color" %in% colnames(bins_pt5)) {
        colors <- unique(bins_pt5$Color)
        names(colors) <- unique(bins_pt5$Var2)
        scale_color <- scale_color_manual(values = colors)
      } else {
        scale_color <- scale_color_viridis(discrete = TRUE)
      }
      
      # plot
      df_melt <- melt(df_log_normalized, variable.factor=TRUE, variable.name="variable", value.name="value")
      bins_pt5 <- table(cut(df_melt$value, breaks=seq(0, 16.5, by=input$bin_size)), df_melt$variable)
      bins_pt5 <- as.data.frame(bins_pt5)
      p <- bins_pt5 %>%
        ggplot(aes(x=Var1, y=Freq, group=Var2, color=Var2, text = paste(ifelse(input$grouped, "Group:", "Sample:"), Var2, "<br>Interval:", Var1, "<br>Count:", Freq))) +
        geom_line() +
        scale_color +
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
      
        ggplotly(p, tooltip = "text")
    }
    
    create_expressed_gene_plot <- function(df, groups, plot_title) {
      numeric_data <- df[,-c(1,2)]
      
      # Function to filter a column and return its length
      filter_column_length <- function(column) {
        # Count only values greater than 5
        return(length(column[column > 0]))
      }
      
      # Get the lengths of filtered columns in a named list
      column_lengths <- sapply(numeric_data, filter_column_length)
      
      # Convert the named list to a data frame with a row for each sample
      lengths_df <- data.frame(Samples = names(column_lengths), Gene_Count = column_lengths)
      
      if (input$grouped) {
        lengths_df_grouped <- lengths_df %>%
          left_join(groups, by = "Samples") 

        color_values <- unique(lengths_df_grouped$Color)
        names(color_values) <- unique(lengths_df_grouped$EXPERIMENTAL_GROUP)
        fill_color <- scale_fill_manual(values = color_values)
        
        p <- ggplot(lengths_df_grouped, aes(x = EXPERIMENTAL_GROUP, y = Gene_Count, fill = EXPERIMENTAL_GROUP)) +
          geom_boxplot(outlier.colour = "red", outlier.shape = 1) +
          fill_color + 
          labs(title = plot_title,
               x = "Groups",
               y = "Gene Count") +
          theme_ipsum() +
          theme(
            axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
            axis.title.x = element_text(vjust = 0.5, hjust = 0.5),
            axis.title.y = element_text(hjust = 0.5),
            panel.border = element_rect(colour = "black", fill = NA, size = 1),
            legend.box.background = element_rect(colour = "black", fill = NA, size = 0.2),
            legend.position = "none",
            legend.title = element_blank(),
            legend.justification = c("right", "top"),
            legend.box.just = "right",
            legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
            legend.key.size = unit(0.5, "cm"))
      } else {
        # Create a bar plot
        p <- ggplot(lengths_df, aes(x = Samples, y = Gene_Count, fill = Samples, text = paste("Sample:", Samples, "<br>Count:", Gene_Count))) +
          geom_bar(stat = "identity") +
          scale_fill_viridis(discrete = TRUE) +
          labs(title = plot_title,
               x = "Samples",
               y = "Gene Count") +
          theme_minimal() +
          theme_ipsum() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_text(vjust = 0.5, hjust = 0.5),
            axis.title.y = element_text(hjust = 0.5),
            panel.border = element_rect(colour = "black", fill = NA, size = 1),
            legend.position = "none",
            legend.box.background = element_rect(colour = "black", fill = NA, size = 0.2),
            legend.title = element_blank(),
            legend.justification = c("right", "top"),
            legend.box.just = "right",
            legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
            legend.key.size = unit(0.5, "cm"))
      }
      
        ggplotly(p, tooltip = "text")
    }
    
    output$hist <- renderPlotly({
      req(dataset$expression_data())
      create_hist_plot()
    })
    
    output$filtered_hist <- renderPlotly({
      req(dataset$expression_data())
      create_filtered_hist_plot()
    })
    
    output$grouped_hist <- renderPlotly({
      req(dataset$expression_data(), dataset$selected_samples_data())
      df <- dataset$expression_data()
      selected_samples <- dataset$selected_samples_data()
      groups <- dataset$groups_data()
      df <- df %>%
        select(c('Symbols', 'Genes', all_of(selected_samples)))
      
      colnames(df)[1] <- "Symbols"
      colnames(df)[2] <- "Genes"
      colnames(groups)[1] <- "Samples"
      colnames(groups)[2] <- "EXPERIMENTAL_GROUP"
      
      create_expressed_gene_plot(df, groups, "Number of Genes before filtering")
    })
    
    
    output$grouped_filtered_hist <- renderPlotly({
      req(dataset$filtered_data())
      df <- dataset$filtered_data()
      groups <- dataset$groups_data()
      
      colnames(df)[1] <- "Symbols"
      colnames(df)[2] <- "Genes"
      colnames(groups)[1] <- "Samples"
      colnames(groups)[2] <- "EXPERIMENTAL_GROUP"
      
      create_expressed_gene_plot(df, groups, "Number of Genes after filtering")
    })
  
    
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
        #filePath4 <- file.path(tempDir, "plot4.png")
        # Save plots
        ggsave(filePath1, plot = create_hist_plot(), device = "png", width = 10, height = 8, dpi = 300)
        ggsave(filePath2, plot = create_filtered_hist_plot(), device = "png", width = 10, height = 8, dpi = 300)
        ggsave(filePath3, plot = create_grouped_plot(), device = "png", width = 10, height = 8, dpi = 300)
        #ggsave(filePath4, plot = create_grouped_filtered_plot(), device = "png", width = 10, height = 8, dpi = 300)
        # Create zip file
        zip(file, files = c(filePath1, filePath2, filePath3, filePath4))
      }
    )
    
  })
}

qcTabUI <- function(id) {
  fluidPage(
    tags$head(
      tags$style(HTML("
      .bottom-centered {
        display: flex;
        align-items: center;
      }
    "))
    ),
    fluidRow(class = "bottom-centered",
      column(2, 
        sliderInput(NS(id, "bin_size"), "Bin size", min = 0, max = 1, value = 0.25, step = 0.25),
        # downloadButton(NS(id, "download_hist"), "Download Plots")
      ),
      column(1, 
        materialSwitch(inputId = NS(id, "grouped"), label = "Show by group: ", value = FALSE, status = "primary")
        # colourInput(NS(id, "color"), "a01_brain_microglia1_Exp1_m:", value = "#FFFFFF"),
        # switchInput(inputId = "grouped", label = "Show by group", value = TRUE)
        # materialSwitch(inputId = "Id077", label = "Show by group", value = TRUE, status = "primary")
      ),
    ),
    fluidRow(
      column(6, 
        plotlyOutput(NS(id, "hist"))
      ),
      column(6, 
        plotlyOutput(NS(id, "filtered_hist"))
      )
    ),
    fluidRow(
      column(6, 
        plotlyOutput(NS(id, "grouped_hist"))
      ),
      column(6, 
        plotlyOutput(NS(id, "grouped_filtered_hist"))
      )
    )
  )
}