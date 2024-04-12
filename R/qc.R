qcTabServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    data <- reactive({
      req(input$expressions)
      
      ext <- tools::file_ext(input$expressions$name)
      switch(ext,
             csv = vroom::vroom(input$expressions$datapath, delim = ","),
             tsv = vroom::vroom(input$expressions$datapath, delim = "\t"),
             validate("Invalid file; Please upload a .csv or .tsv file")
      )
    
    })
    
    groups <- reactive({
      req(input$groups)
      
      ext <- tools::file_ext(input$groups$name)
      switch(ext,
             csv = vroom::vroom(input$groups$datapath, delim = ","),
             tsv = vroom::vroom(input$groups$datapath, delim = "\t"),
             validate("Invalid file; Please upload a .csv or .tsv file")
      )
    })
    
    grouped_data <- reactive({
      req(input$expressions)
      req(input$groups)
      
      df <- data.frame(data())
      
      new_col_name <- c(colnames(df)[1], colnames(df)[2])
      
      for (old in colnames(df)[3:ncol(df)]) {
        new_col_name <- append(new_col_name, paste0(groups()[which(groups()$SAMPLES == gsub("\\.", " ", old)), 2]))
      }
    
      colnames(df) <- new_col_name
      
      # 找出所有重复的列名的索引
      dup <- names(df) %in% names(df)[duplicated(names(df))]
      indices <- which(dup)
      
      # 按照同名的列分组
      groups <- split(indices, names(df)[indices])
      
      # 对每一组同名列进行相加
      new_cols <- sapply(groups, function(idxs) rowSums(df[, idxs]))
      
      # 找出所有不重复的列
      non_dup_cols <- df[, !duplicated(names(df))]
      
      # 将结果合并为一个新的数据框
      new_df <- cbind(non_dup_cols, as.data.frame(new_cols))
      
      df
    })
    
    output$group_table<- renderTable({
      #browser()
      head(grouped_data(), 10)
    })
    
    
    output$head <- renderTable({
      #browser()
      head(data(), 10)
    })
    
    output$hist <- renderPlot({
      #browser()
      df <- data()
      df_log_normalized <- log1p(df[3:ncol(df)])
      df_vector <- unlist(df_log_normalized)
      #browser()
      bins <- seq(from = min(df_vector), to = max(df_vector), length.out = max(df_vector) / input$bin_size)
      hist(df_vector, breaks = bins, xlab = paste("Bins(binsize=", input$bin_size, ")"), ylab = "Gene Counts", main = "")
      #ggplot(data()[,2], mapping = aes(x = "a01 brain microglia1 Exp1 m", y = "a01 brain microglia1 Exp1 m"))
    }, res = 96)
    
    output$grouped_hist <- renderPlot({
      #browser()
      df <- grouped_data()
      df <- data()
      df_log_normalized <- log1p(df[3:ncol(df)])
      df_vector <- unlist(df_log_normalized)
      #browser()
      bins <- seq(from = min(df_log_normalized), to = max(df_log_normalized), length.out = max(df_log_normalized) / input$bin_size)
      hist(df_vector, breaks = bins, xlab = paste("Bins(binsize=", input$bin_size, ")"), ylab = "Gene Counts", main = "")
      #ggplot(data()[,2], mapping = aes(x = "a01 brain microglia1 Exp1 m", y = "a01 brain microglia1 Exp1 m"))
    }, res = 96)
    
    observeEvent(data(), {
      x <- log2(data()$`b05 Neutrophil1 Exp2` + 1)
      updateSliderInput(inputId = "bin_size", min = min(x), max = max(x) %/% 2)
    })  
  })
}

qcTabUI <- function(id) {
  fluidPage(
    titlePanel("QC"),
    fluidRow(
      column(6, 
        fileInput(NS(id, "expressions"), NULL, accept = c(".csv", ".tsv")),
      ),
      column(6, 
        fileInput(NS(id, "groups"), NULL, accept = c(".csv", ".tsv")),
      )
    ),
    fluidRow(
      column(6, 
        sliderInput(NS(id, "bin_size"), "Bin size", min = 0, max = 3, value = 1, step = 0.01),
        plotOutput(NS(id, "hist"))
      ),
      column(6, 
        sliderInput(NS(id, "bin_size"), "Bin size", min = 0, max = 3, value = 1, step = 0.01),
        plotOutput(NS(id, "grouped_hist"))
      )
    )
  )
}