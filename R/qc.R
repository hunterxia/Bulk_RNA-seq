library(viridis)
library(hrbrthemes)
library(ggplot2)
library(tidyr)
library(shinyWidgets)
library(plotly)
library(logger)
library(data.table)
library(future)
library(logger)

qcTabServer <- function(id, dataset) {
  moduleServer(id, function(input, output, session) {

    expression_data_cache <- reactiveVal(NULL)

    observe({
      log_info("Loading expression data into cache asynchronously")
      expression_data_cache(dataset$expression_data())
    })

    grouped <- reactive({ TRUE })

    get_max_freq <- function(df, groups) {
      #df_long <- pivot_longer(df, cols = -c(Symbol, Gene_Symbol), names_to = "Samples", values_to = "Values")
      df_long <- as.data.table(df) %>%
        melt(id.vars = c("Symbol", "Gene_Symbol"), variable.name = "Samples", value.name = "Values")
      df_long <- df_long %>%
        left_join(groups, by = "Samples")
      df_long$LogValue <- log1p(df_long$Values)
      # df_long$Bins <- cut(df_long$LogValue, breaks = seq(0, 16.5, by = input$bin_size), right = FALSE)

      df_long$Bins <- cut(df_long$LogValue,
                          breaks = seq(0, max(df_long$LogValue, na.rm = TRUE),
                                       by = input$bin_size),
                          right = TRUE,
                          include.lowest = FALSE)

      bins_df <- df_long %>%
        group_by(Bins, Samples, EXPERIMENTAL_GROUP) %>%
        summarise(Freq = n(), .groups = 'drop')

      return(max(bins_df$Freq))
    }

    create_hist_plot <- function() {
      # pdf(NULL)

      df <- expression_data_cache()
      groups <- dataset$groups_data()
      selected_samples <- dataset$selected_samples_data()

      df <- df %>%
        select(c('Symbol', 'Gene_Symbol', all_of(selected_samples)))

      colnames(df)[1] <- "Symbol"
      colnames(df)[2] <- "Gene_Symbol"
      colnames(groups)[1] <- "Samples"
      colnames(groups)[2] <- "EXPERIMENTAL_GROUP"

      # Pivot data to long format
      #df_long <- pivot_longer(df, cols = -c(Symbol, Gene_Symbol), names_to = "Samples", values_to = "Values")
      df_long <- as.data.table(df) %>%
        melt(id.vars = c("Symbol", "Gene_Symbol"), variable.name = "Samples", value.name = "Values")

      # Join groups data to assign experimental groups
      df_long <- df_long %>%
        left_join(groups, by = "Samples")

      # remove all rows with count is 0
      df_long <- df_long[df_long$Values != 0,]

      # Log normalize the values
      df_long$LogValue <- log1p(df_long$Values)

      # Cut the LogValue into bins
      df_long$Bins <- cut(df_long$LogValue,
                          breaks = seq(0, max(df_long$LogValue, na.rm = TRUE),
                                       by = input$bin_size),
                          right = TRUE,
                          include.lowest = FALSE)


      # Calculate frequencies
      bins_df <- df_long %>%
        group_by(Bins, Samples, EXPERIMENTAL_GROUP) %>%
        summarise(Freq = n(), .groups = 'drop')

      # Assign colors based on input$grouped
      # if (TRUE) {
      #   color_mapping <- bins_df$EXPERIMENTAL_GROUP
      #   color_label <- "Experimental Group"
      #
      #   if ("Color" %in% colnames(groups)) {
      #     color_values <- groups %>%
      #       distinct(EXPERIMENTAL_GROUP, Color) %>%
      #       deframe()
      #     scale_color <- scale_color_manual(values = color_values)
      #   } else {
      #     scale_color <- scale_color_viridis(discrete = TRUE)
      #   }
      #
      # } else {
      #   color_mapping <- bins_df$Samples
      #   color_label <- "Samples"
      #   scale_color <- scale_color_viridis(discrete = TRUE)
      # }
      # #max_freq <- get_max_freq(dataset$expression_data(), groups)
      # max_freq <- get_max_freq(expression_data_cache(), groups)
      #
      # # Plot
      # p <- bins_df %>%
      #   ggplot(aes(x = Bins, y = Freq, group = Samples, color = color_mapping,
      #              text = paste("Sample:", Samples, "<br>Group:", EXPERIMENTAL_GROUP, "<br>Interval:", Bins, "<br>Count:", Freq))) +
      #   geom_line() +
      #   scale_color +
      #   ggtitle("Histogram before filtering") +
      #   theme_ipsum() +
      #   labs(
      #     x = NULL,
      #     y = "Gene Counts",
      #     color = NULL
      #   ) +
      #   #scale_y_continuous(limits = c(0, max_freq))+
      #   scale_y_continuous(limits = c(0, 4000)) +
      #   guides(color = guide_legend(nrow = 6, byrow = TRUE, title.position = "top")) +
      #   theme(
      #     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      #     axis.title.x = element_text(vjust = 0.5, hjust = 0.5),
      #     axis.title.y = element_text(hjust = 0.5),
      #     panel.border = element_rect(colour = "black", fill = NA, size = 1),
      #     legend.box.background = element_rect(colour = "black", fill = NA, size = 0.2),
      #     legend.position = c(1, 1),
      #     legend.justification = c("right", "top"),
      #     legend.box.just = "right",
      #     legend.margin = margin(0, 0, 0, 0, "pt"),
      #     legend.key.size = unit(0.5, "cm")
      #   )
      # # ggplotly(p, tooltip = "text") %>%
      # #   layout(
      # #     width = 400,
      # #     height = 400
      # #   )
      # plotly_obj <- ggplotly(p, tooltip = "text") %>%
      #   layout(
      #     width = 400,
      #     height = 400
      #   )
      # dev.off()

      custom_colors <- NULL
      if ("Color" %in% colnames(groups)) {
        color_values <- groups %>%
          distinct(EXPERIMENTAL_GROUP, Color) %>%
          deframe()
        custom_colors <- unname(color_values)
      }

      plotly_obj <- plot_ly(
        data = bins_df,
        x = ~Bins,
        y = ~Freq,
        type = 'scatter',
        mode = 'lines',
        split = ~Samples,
        color = ~EXPERIMENTAL_GROUP,
        colors = if (!is.null(custom_colors)) custom_colors else "viridis",
        text = ~paste("Sample:", Samples,
                      "<br>Group:", EXPERIMENTAL_GROUP,
                      "<br>Interval:", Bins,
                      "<br>Count:", Freq),
        hoverinfo = "text"
      ) %>%
        layout(
          title = "Histogram before filtering",
          xaxis = list(title = NULL, tickangle = 90),
          yaxis = list(title = "Gene Counts", range = c(0, 4000)),
          # width = 400,
          # height = 400,
          legend = list(orientation = "v", x = 1, y = 1)
        )
      plotly_obj
    }

    create_filtered_hist_plot <- function() {
      # pdf(NULL)

      df <- dataset$filtered_data()
      groups <- dataset$groups_data()

      colnames(df)[1] <- "Symbol"
      colnames(df)[2] <- "Gene_Symbol"
      colnames(groups)[1] <- "Samples"
      colnames(groups)[2] <- "EXPERIMENTAL_GROUP"

      # Pivot data to long format
      #df_long <- pivot_longer(df, cols = -c(Symbol, Gene_Symbol), names_to = "Samples", values_to = "Values")
      df_long <- as.data.table(df) %>%
        melt(id.vars = c("Symbol", "Gene_Symbol"), variable.name = "Samples", value.name = "Values")

      # Join groups data to assign experimental groups
      df_long <- df_long %>%
        left_join(groups, by = "Samples")

      # remove all rows with count is 0
      df_long <- df_long[df_long$Values != 0,]

      # Log normalize the values
      df_long$LogValue <- log1p(df_long$Values)

      # Cut the LogValue into bins
      df_long$Bins <- cut(df_long$LogValue,
                          breaks = seq(0, max(df_long$LogValue, na.rm = TRUE),
                                       by = input$bin_size),
                          right = TRUE,
                          include.lowest = FALSE)


      # Calculate frequencies
      bins_df <- df_long %>%
        group_by(Bins, Samples, EXPERIMENTAL_GROUP) %>%
        summarise(Freq = n(), .groups = 'drop')

      # # Assign colors based on input$grouped
      # if (TRUE) {
      #   color_mapping <- bins_df$EXPERIMENTAL_GROUP
      #   color_label <- "Experimental Group"
      #
      #   if ("Color" %in% colnames(groups)) {
      #     color_values <- groups %>%
      #       distinct(EXPERIMENTAL_GROUP, Color) %>%
      #       deframe()
      #     scale_color <- scale_color_manual(values = color_values)
      #   } else {
      #     scale_color <- scale_color_viridis(discrete = TRUE)
      #   }
      #
      # } else {
      #   color_mapping <- bins_df$Samples
      #   color_label <- "Samples"
      #   scale_color <- scale_color_viridis(discrete = TRUE)
      # }

      #max_freq <- get_max_freq(dataset$expression_data(), groups)
      max_freq <- get_max_freq(expression_data_cache(), groups)

      # Plot
      # p <- bins_df %>%
      #   ggplot(aes(x = Bins, y = Freq, group = Samples, color = color_mapping,
      #              text = paste("Sample:", Samples, "<br>Group:", EXPERIMENTAL_GROUP, "<br>Interval:", Bins, "<br>Count:", Freq))) +
      #   geom_line() +
      #   scale_color +
      #   ggtitle("Histogram after filtering") +
      #   theme_ipsum() +
      #   labs(
      #     x = NULL,
      #     y = "Gene Counts",
      #     color = NULL
      #   ) +
      #   #scale_y_continuous(limits = c(0, max_freq)) +
      #   scale_y_continuous(limits = c(0, 4000)) +
      #   guides(color = guide_legend(nrow = 6, byrow = TRUE, title.position = "top")) +
      #   theme(
      #     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      #     axis.title.x = element_text(vjust = 0.5, hjust = 0.5),
      #     axis.title.y = element_text(hjust = 0.5),
      #     panel.border = element_rect(colour = "black", fill = NA, size = 1),
      #     legend.box.background = element_rect(colour = "black", fill = NA, size = 0.2),
      #     legend.position = c(1, 1),
      #     legend.justification = c("right", "top"),
      #     legend.box.just = "right",
      #     legend.margin = margin(0, 0, 0, 0, "pt"),
      #     legend.key.size = unit(0.5, "cm")
      #   )
      # # ggplotly(p, tooltip = "text") %>%
      # #   layout(
      # #     width = 400,
      # #     height = 400
      # #   )
      # plotly_obj <- ggplotly(p, tooltip = "text") %>%
      #   layout(
      #     width = 400,
      #     height = 400
      #   )

      custom_colors <- NULL
      if ("Color" %in% colnames(groups)) {
        color_values <- groups %>%
          distinct(EXPERIMENTAL_GROUP, Color) %>%
          deframe()
        custom_colors <- unname(color_values)
      }

      plotly_obj <- plot_ly(
        data = bins_df,
        x = ~Bins,
        y = ~Freq,
        type = 'scatter',
        mode = 'lines',
        split = ~Samples,
        color = ~EXPERIMENTAL_GROUP,
        colors = if (!is.null(custom_colors)) custom_colors else "viridis",
        text = ~paste("Sample:", Samples,
                      "<br>Group:", EXPERIMENTAL_GROUP,
                      "<br>Interval:", Bins,
                      "<br>Count:", Freq),
        hoverinfo = "text"
      ) %>%
        layout(
          title = "Histogram after filtering",
          xaxis = list(title = ""),
          yaxis = list(title = "Gene Counts", range = c(0, 4000)),
          # width = 400,
          # height = 400,
          legend = list(orientation = "h", x = 1, y = 1)
        )

      # dev.off()
      plotly_obj
    }

    create_expressed_gene_plot <- function(df, groups, plot_title) {
      # pdf(NULL)

      numeric_data <- df[, -c(1, 2)]

      # Function to count non-zero values per sample
      count_nonzero <- function(column) {
        sum(column > 0)
      }

      # Get the counts of non-zero values per sample
      counts <- sapply(numeric_data, count_nonzero)

      counts_df <- data.frame(Samples = names(counts), Gene_Count = counts)

      # Join groups data to assign experimental groups
      counts_df <- counts_df %>%
        left_join(groups, by = "Samples")

      # if (TRUE) {
      #   color_mapping <- counts_df$EXPERIMENTAL_GROUP
      #   color_label <- "Experimental Group"
      #
      #   if ("Color" %in% colnames(groups)) {
      #     color_values <- groups %>%
      #       distinct(EXPERIMENTAL_GROUP, Color) %>%
      #       deframe()
      #     fill_color <- scale_fill_manual(values = color_values)
      #   } else {
      #     fill_color <- scale_fill_viridis(discrete = TRUE)
      #   }
      #
      # } else {
      #   color_mapping <- counts_df$Samples
      #   color_label <- "Samples"
      #   fill_color <- scale_fill_viridis(discrete = TRUE)
      # }
      # #unfiltered_max <- max(sapply(dataset$expression_data()[, -c(1, 2)], function(column) sum(column > 0)))
      # unfiltered_max <- max(sapply(expression_data_cache()[, -c(1, 2)], function(column) sum(column > 0)))
      #
      # # Plot
      # p <- ggplot(counts_df, aes(x = Samples, y = Gene_Count, fill = color_mapping,
      #                            text = paste("Sample:", Samples, "<br>Group:", EXPERIMENTAL_GROUP, "<br>Count:", Gene_Count))) +
      #   geom_bar(stat = "identity") +
      #   fill_color +
      #   labs(title = plot_title,
      #        x = NULL,
      #        y = "Gene Count",
      #        fill = NULL) +
      #   theme_ipsum() +
      #   scale_y_continuous(limits = c(0, unfiltered_max)) +
      #   theme(
      #     axis.text.x = element_text(angle = 45, hjust = 1),
      #     axis.title.x = element_text(vjust = 0.5, hjust = 0.5),
      #     axis.title.y = element_text(hjust = 0.5),
      #     panel.border = element_rect(colour = "black", fill = NA, size = 1),
      #     legend.position = "right",
      #     legend.box.background = element_rect(colour = "black", fill = NA, size = 0.2),
      #     legend.title = element_blank(),
      #     legend.justification = c("right", "top"),
      #     legend.box.just = "right",
      #     legend.margin = margin(0, 0, 0, 0, "pt"),
      #     legend.key.size = unit(0.5, "cm")
      #   )
      # # ggplotly(p, tooltip = "text") %>%
      # #   layout(
      # #     width = 400,
      # #     height = 400
      # #   )
      # plotly_obj <- ggplotly(p, tooltip = "text") %>%
      #   layout(
      #     width = 400,
      #     height = 400
      #   )

      custom_colors <- NULL
      if ("Color" %in% colnames(groups)) {
        color_values <- groups %>%
          distinct(EXPERIMENTAL_GROUP, Color) %>%
          deframe()
        custom_colors <- unname(color_values)
      }

      unfiltered_max <- max(sapply(expression_data_cache()[, -c(1, 2)], function(column) sum(column > 0)))

      plotly_obj <- plot_ly(
        data = counts_df,
        x = ~Samples,
        y = ~Gene_Count,
        type = 'bar',
        text = ~paste("Sample:", Samples,
                      "<br>Group:", EXPERIMENTAL_GROUP,
                      "<br>Count:", Gene_Count),
        hoverinfo = "text",
        color = ~EXPERIMENTAL_GROUP,
        colors = if (!is.null(custom_colors)) custom_colors else "viridis"
      ) %>%
        layout(
          title = plot_title,
          xaxis = list(title = "", tickangle = -45),
          yaxis = list(title = "Gene Count", range = c(0, unfiltered_max)),
          # width = 400,
          # height = 400,
          legend = list(orientation = "v")
        )
      # dev.off()
      plotly_obj
    }

    output$hist <- renderPlotly({
      #req(dataset$expression_data())
      req(expression_data_cache())
      create_hist_plot()
    })

    output$filtered_hist <- renderPlotly({
      #req(dataset$expression_data())
      req(expression_data_cache())
      create_filtered_hist_plot()
    })

    output$grouped_hist <- renderPlotly({
      req(dataset$expression_data(), dataset$selected_samples_data())
      #df <- dataset$expression_data()
      df <- expression_data_cache()
      selected_samples <- dataset$selected_samples_data()
      groups <- dataset$groups_data()
      df <- df %>%
        select(c('Symbol', 'Gene_Symbol', all_of(selected_samples)))

      colnames(df)[1] <- "Symbol"
      colnames(df)[2] <- "Gene_Symbol"
      colnames(groups)[1] <- "Samples"
      colnames(groups)[2] <- "EXPERIMENTAL_GROUP"

      create_expressed_gene_plot(df, groups, "Number of Genes before filtering")
    })


    output$grouped_filtered_hist <- renderPlotly({
      req(dataset$filtered_data())
      df <- dataset$filtered_data()
      groups <- dataset$groups_data()

      colnames(df)[1] <- "Symbol"
      colnames(df)[2] <- "Gene_Symbol"
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
      tags$style(HTML(" .bottom-centered { display: flex; align-items: center; } "))
    ),
    fluidRow(class = "bottom-centered",
             column(2,
                    sliderInput(NS(id, "bin_size"), "Bin size", min = 0, max = 1, value = 0.25, step = 0.05)
             )
    ),
    fluidRow(
      column(6,
             plotlyOutput(NS(id, "hist"), width = "100%", height = "100%")
      ),
      column(6,
             plotlyOutput(NS(id, "filtered_hist"), width = "100%", height = "100%")
      )
    ),
    fluidRow(
      column(6,
             plotlyOutput(NS(id, "grouped_hist"), width = "100%", height = "100%")
      ),
      column(6,
             plotlyOutput(NS(id, "grouped_filtered_hist"), width = "100%", height = "100%")
      )
    )
  )
}