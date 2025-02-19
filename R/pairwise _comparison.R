library(shiny)
library(ggplot2)
library(plotly)
library(DT)
library(dplyr)
library(DESeq2)
library(logger)

PairwiseComparisonTabUI <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    titlePanel("Pairwise Comparison of Gene Expression"),
    sidebarLayout(
      sidebarPanel(
        selectInput(ns("select_group1"), "Select Group 1", choices = NULL),
        selectInput(ns("select_group2"), "Select Group 2", choices = NULL),
        selectInput(ns("test_type"), "Select Test Type",
                    choices = c(
                      "Unpaired t-test" = "paired",  
                      #"Unpaired t-test" = "unpaired",
                      "Deseq" = "deseq"
                    )
        ),
        numericInput(
          ns("fc_cutoff"),
          "Fold Change Cutoff",
          value = 2,
          step = 0.05
        ),
        numericInput(
          ns("pvalue_cutoff"),
          "P-value Cutoff",
          value = 0.5,
          min = 0,
          max = 1,
          step = 0.01
        ),
        checkboxInput(ns("adjust_pvalue"), "Use Adjusted P-Value", value = FALSE),
        actionButton(ns("run_analysis"), "Run Analysis"),
        downloadButton(ns("download_all"), "Download Full DEG Output"),
        downloadButton(ns("download_filtered"), "Download Filtered DEGs")
      ),
      mainPanel(
        tabsetPanel(type = "tabs",
                    tabPanel("Volcano Plot",
                             textOutput(ns("volcano_gene_counts")),
                             plotlyOutput(ns("volcano_plot"), width = "100%", height = "100%")),
                    tabPanel("MA Plot",
                             textOutput(ns("ma_gene_counts")),
                             plotlyOutput(ns("ma_plot"), width = "100%", height = "100%")),
                    tabPanel("Positive DEG Table", DTOutput(ns("pos_deg_table"))),
                    tabPanel("Negative DEG Table", DTOutput(ns("neg_deg_table")))
        )
      )
    )
  )
}

perform_DEG_analysis <- function(exp_data, group1_samples, group2_samples, test_type, fc_cutoff, pvalue_cutoff, adjust_pvalue) {
  log_debug("Starting DEG analysis")
  log_trace("Input parameters: FC cutoff={fc_cutoff}, P-value={pvalue_cutoff}")
  
  if (test_type == "deseq") {
    exp_data_filtered <- exp_data[, c(group1_samples, group2_samples)]
    exp_data_filtered <- round(exp_data_filtered)
    
    condition <- factor(c(rep("group1", length(group1_samples)), rep("group2", length(group2_samples))))
    
    tryCatch({
      dds <- DESeqDataSetFromMatrix(
        countData = exp_data_filtered,
        colData = data.frame(condition = condition),
        design = ~condition
      )
      log_debug("DESeq object created with {ncol(dds)} samples")
      
      dds <- DESeq(dds)
      log_debug("DESeq analysis completed: {length(resultsNames(dds))} coefficients")
      
      res <- results(dds)
      log_trace("DESeq results range: log2FC [{min(res$log2FoldChange)}, {max(res$log2FoldChange)}]")
    }, error = function(e) {
      log_error("DESeq2 failed: {e$message}")
      stop(e)
    })
    
    log2_fold_changes <- res$log2FoldChange
    p_values <- res$pvalue
    adj_pvalues <- res$padj  # Always get adjusted p-values from DESeq2
    
    deg_results <- data.frame(
      Symbols = rownames(exp_data),
      Gene_Symbol = exp_data$Gene_Symbol,
      baseMean = res$baseMean,
      Group1Mean = rowMeans(exp_data_filtered[, group1_samples], na.rm = TRUE),
      Group2Mean = rowMeans(exp_data_filtered[, group2_samples], na.rm = TRUE),
      Log2FoldChange = log2_fold_changes,
      lfcSE = res$lfcSE,
      stat = res$stat,
      PValue = p_values,
      AdjPValue = adj_pvalues,
      stringsAsFactors = FALSE
    )
    
  } else {
    # Filter expression data to include only the relevant samples
    group1_data <- exp_data[, colnames(exp_data) %in% group1_samples, drop = FALSE]
    group2_data <- exp_data[, colnames(exp_data) %in% group2_samples, drop = FALSE]
    
    # Add pseudocount to avoid log transformation issues
    pseudocount <- 1
    group1_data <- group1_data + pseudocount
    group2_data <- group2_data + pseudocount
    
    # Ensure the data is numeric
    group1_data <- data.frame(lapply(group1_data, as.numeric))
    group2_data <- data.frame(lapply(group2_data, as.numeric))
    
    # Calculate log2 fold changes
    log2_fold_changes <- log2(rowMeans(group2_data, na.rm = TRUE)) - log2(rowMeans(group1_data, na.rm = TRUE))
    
    # Calculate p-values for each gene
    p_values <- vector("numeric", length = nrow(exp_data))
    for (i in 1:nrow(exp_data)) {
      gene_data1 <- as.numeric(group1_data[i, , drop = TRUE])
      gene_data2 <- as.numeric(group2_data[i, , drop = TRUE])
      
      if (test_type == "paired") {
         p_values[i] <- t.test(gene_data1, gene_data2, paired = TRUE)$p.value
       } else if (test_type == "unpaired") {
        p_values[i] <- t.test(gene_data1, gene_data2, paired = FALSE)$p.value
      } else if (test_type == "mann_whitney") {
        p_values[i] <- wilcox.test(gene_data1, gene_data2, exact = FALSE)$p.value
      } else if (test_type == "wilcoxon") {
        p_values[i] <- wilcox.test(gene_data1, gene_data2, paired = TRUE, exact = FALSE)$p.value
      }
    }
    
    # Always adjust p-values
    adj_pvalues <- p.adjust(p_values, method = "BH")
    
    deg_results <- data.frame(
      Symbols = rownames(exp_data),
      Gene_Symbol = exp_data$Gene_Symbol,
      Group1Mean = rowMeans(group1_data, na.rm = TRUE),
      Group2Mean = rowMeans(group2_data, na.rm = TRUE),
      Log2FoldChange = log2_fold_changes,
      PValue = p_values,
      AdjPValue = adj_pvalues,
      stringsAsFactors = FALSE
    )
  }
  
  # Determine which p-values to use for significance
  if (adjust_pvalue) {
    pvalues_for_significance <- deg_results$AdjPValue
  } else {
    pvalues_for_significance <- deg_results$PValue
  }
  
  # Determine significance based on the user-defined cutoffs
  deg_results$Significant <- (pvalues_for_significance <= pvalue_cutoff) &
    (abs(deg_results$Log2FoldChange) >= log2(fc_cutoff))
  
  # Categorize regulation status
  deg_results$Regulation <- "Not Significant"
  deg_results$Regulation[deg_results$Significant & deg_results$Log2FoldChange >= log2(fc_cutoff)] <- "Upregulated"
  deg_results$Regulation[deg_results$Significant & deg_results$Log2FoldChange <= -log2(fc_cutoff)] <- "Downregulated"
  
  return(deg_results)
}

generate_volcano_plot <- function(results, group1_name, group2_name) {
  total_genes <- nrow(results)
  up_genes <- sum(results$Regulation == "Upregulated", na.rm = TRUE)
  down_genes <- sum(results$Regulation == "Downregulated", na.rm = TRUE)
  
  plotly_obj <- plot_ly(
    data = results,
    x = ~Log2FoldChange,
    y = ~-log10(PValue),
    type = 'scatter',
    mode = 'markers',
    color = ~Regulation,
    colors = c("Upregulated" = "red",
               "Downregulated" = "blue",
               "Not Significant" = "grey"),
    text = ~paste(
      "Gene:", Gene_Symbol,
      "<br>", group1_name, "Mean:", Group1Mean,
      "<br>", group2_name, "Mean:", Group2Mean,
      "<br>Raw Fold Change:", round(2^Log2FoldChange, 2),
      "<br>Log2 Fold Change:", round(Log2FoldChange, 2),
      "<br>P-value:", signif(PValue, 3),
      "<br>Adj P-value:", signif(AdjPValue, 3)
    ),
    hoverinfo = "text",
    marker = list(opacity = 0.8)
  ) %>%
    layout(
      title = "Volcano Plot",
      xaxis = list(title = "Log2 Fold Change"),
      yaxis = list(title = "-log10(P-value)")
    )
  
  plotly_obj
}

generate_ma_plot <- function(exp_data, group1_samples, group2_samples, results, group1_name, group2_name) {
  pseudocount <- 1
  group1_data <- exp_data[, colnames(exp_data) %in% group1_samples, drop = FALSE] + pseudocount
  group2_data <- exp_data[, colnames(exp_data) %in% group2_samples, drop = FALSE] + pseudocount
  
  A <- (log2(rowMeans(group1_data)) + log2(rowMeans(group2_data))) / 2
  M <- log2(rowMeans(group2_data)) - log2(rowMeans(group1_data))
  
  ma_data <- data.frame(
    A = A,
    M = M,
    Regulation = results$Regulation,
    Symbols = results$Symbols,
    Gene_Symbol = results$Gene_Symbol,
    Group1Mean = rowMeans(group1_data),
    Group2Mean = rowMeans(group2_data),
    Log2FoldChange = M,
    PValue = results$PValue,
    AdjPValue = results$AdjPValue,
    stringsAsFactors = FALSE
  )
  
  total_genes <- nrow(ma_data)
  up_genes <- sum(ma_data$Regulation == "Upregulated", na.rm = TRUE)
  down_genes <- sum(ma_data$Regulation == "Downregulated", na.rm = TRUE)
  
  print(paste("Total genes in MA plot:", total_genes))
  print(paste("Upregulated genes:", up_genes))
  print(paste("Downregulated genes:", down_genes))
  
  plotly_obj <- plot_ly(
    data = ma_data,
    x = ~A,
    y = ~M,
    type = 'scatter',
    mode = 'markers',
    color = ~Regulation,
    colors = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey"),
    text = ~paste(
      "Gene:", Gene_Symbol,
      "<br>", group1_name, "Mean:", Group1Mean,
      "<br>", group2_name, "Mean:", Group2Mean,
      "<br>Raw Fold Change:", round(2^Log2FoldChange, 2),
      "<br>Log2 Fold Change:", round(Log2FoldChange, 2),
      "<br>P-value:", signif(PValue, 3),
      "<br>Adj P-value:", signif(AdjPValue, 3)
    ),
    hoverinfo = "text",
    marker = list(opacity = 0.5)
  ) %>%
    layout(
      title = "MA Plot",
      xaxis = list(title = "Average Log Intensity (A)"),
      yaxis = list(title = "Log Fold Change (M)")
    )
  
  plotly_obj
}

PairwiseComparisonTabServer <- function(id, dataset) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    log_info("Initializing Pairwise Comparison Tab")
    
    # Access processed data from the dataset reactive
    expression_data <- reactive(dataset$filtered_data())
    groups_data <- reactive(dataset$groups_data())
    selected_groups <- reactive(dataset$selected_groups())
    
    # Store the original DEG analysis results when analysis is run
    original_deg_results <- reactiveVal()
    
    # Disable download buttons initially
    shinyjs::disable("download_all")
    shinyjs::disable("download_filtered")
    
    observe({
      req(selected_groups())
      log_debug("Updating group selections: {paste(selected_groups(), collapse = ', ')}")
      grp_names <- selected_groups()
      if (is.null(grp_names)) return()
      updateSelectInput(session, "select_group1", choices = grp_names, selected = grp_names[1])
      updateSelectInput(session, "select_group2", choices = grp_names, selected = grp_names[length(grp_names)])
    })
    
    observeEvent(input$run_analysis, {
      req(expression_data(), groups_data())
      exp_data <- expression_data()
      grp_data <- groups_data()
      log_info("Running DEG analysis")
      log_debug("Parameters: Group1={input$select_group1}, Group2={input$select_group2}, Test={input$test_type}")
      start_time <- Sys.time()
      
      showModal(modalDialog("Running analysis, please wait...", footer = NULL))
      
      tryCatch({
        # Map group selections to sample names
        group1_samples <- grp_data$Sample[grp_data$Group == input$select_group1]
        group2_samples <- grp_data$Sample[grp_data$Group == input$select_group2]
        group1_samples <- intersect(group1_samples, colnames(exp_data))
        group2_samples <- intersect(group2_samples, colnames(exp_data))
        
        log_info("Starting DEG analysis ({input$test_type})")
        # Perform DEG analysis
        res <- perform_DEG_analysis(
          exp_data, group1_samples, group2_samples,
          input$test_type, input$fc_cutoff, input$pvalue_cutoff, input$adjust_pvalue
        )
        
        # Store the original results for later re-evaluation of significance
        original_deg_results(res)
        
        shinyjs::enable("download_all")
        shinyjs::enable("download_filtered")
        
        # Update DEG tables and gene counts
        pos_results <- res[res$Regulation == "Upregulated", ]
        neg_results <- res[res$Regulation == "Downregulated", ]
        total_genes <- nrow(res)
        up_genes <- sum(res$Regulation == "Upregulated", na.rm = TRUE)
        down_genes <- sum(res$Regulation == "Downregulated", na.rm = TRUE)
        
        output$volcano_gene_counts <- renderText({
          paste("Total genes:", total_genes,
                "| Upregulated genes:", up_genes,
                "| Downregulated genes:", down_genes)
        })
        output$ma_gene_counts <- renderText({
          paste("Total genes:", total_genes,
                "| Upregulated genes:", up_genes,
                "| Downregulated genes:", down_genes)
        })
        output$pos_deg_table <- renderDT({
          datatable(pos_results, options = list(pageLength = 10), rownames = FALSE)
        })
        output$neg_deg_table <- renderDT({
          datatable(neg_results, options = list(pageLength = 10), rownames = FALSE)
        })
      }, error = function(e) {
        log_error("Analysis failed: {e$message}")
      }, finally = {
        log_info("Total execution time: {round(Sys.time() - start_time, 2)}s")
        removeModal()
      })
    })
    
    # Create a reactive expression that recalculates significance based on the current cutoffs
    reactive_deg_results <- reactive({
      req(original_deg_results())
      results <- original_deg_results()
      p_val <- if (input$adjust_pvalue) results$AdjPValue else results$PValue
      
      results$Significant <- (p_val <= input$pvalue_cutoff) &
        (abs(results$Log2FoldChange) >= log2(input$fc_cutoff))
      
      results$Regulation <- "Not Significant"
      results$Regulation[results$Significant & results$Log2FoldChange >= log2(input$fc_cutoff)] <- "Upregulated"
      results$Regulation[results$Significant & results$Log2FoldChange <= -log2(input$fc_cutoff)] <- "Downregulated"
      
      results
    })
    
    # Automatically update volcano plot when cutoffs change
    output$volcano_plot <- renderPlotly({
      req(reactive_deg_results())
      generate_volcano_plot(reactive_deg_results(), input$select_group1, input$select_group2)
    })
    
    # Automatically update MA plot when cutoffs change
    output$ma_plot <- renderPlotly({
      req(expression_data(), groups_data(), reactive_deg_results())
      exp_data <- expression_data()
      grp_data <- groups_data()
      group1_samples <- grp_data$Sample[grp_data$Group == input$select_group1]
      group2_samples <- grp_data$Sample[grp_data$Group == input$select_group2]
      group1_samples <- intersect(group1_samples, colnames(exp_data))
      group2_samples <- intersect(group2_samples, colnames(exp_data))
      
      generate_ma_plot(exp_data, group1_samples, group2_samples, reactive_deg_results(), input$select_group1, input$select_group2)
    })
    
    # Download handlers for DEG results
    output$download_all <- downloadHandler(
      filename = function() { "full_deg_output.csv" },
      content = function(file) {
        req(reactive_deg_results())
        write.csv(reactive_deg_results(), file, row.names = FALSE)
      }
    )
    output$download_filtered <- downloadHandler(
      filename = function() { "filtered_deg_output.csv" },
      content = function(file) {
        req(reactive_deg_results())
        write.csv(subset(reactive_deg_results(), Regulation != "Not Significant"), file, row.names = FALSE)
      }
    )
  })
}
