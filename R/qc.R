qcTabServer <- function(id, dataset) {
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
          x = paste0("Bins(binsize=", input$bin_size, ")"),  # X轴名称
          y = "Gene Counts",          # Y轴名称
          title = "Histogram of Gene Expression"
        ) +
        guides(color = guide_legend(nrow = 6, byrow = TRUE, title.position = "top")) +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title.x = element_text(vjust = 0.5, hjust = 0.5),  # X轴标题垂直居中
          axis.title.y = element_text(hjust = 0.5),  # Y轴标题水平居中
          panel.border = element_rect(colour = "black", fill = NA, size = 1),  # 给图表加上边框
          legend.box.background = element_rect(colour = "black", fill = NA, size = 0.2),  # 给图例加上边框
          legend.title = element_blank(),
          legend.position = c(1, 1),  # 将图例放在图表内部，位置根据需要调整
          legend.justification = c("right", "top"),  # 图例的对齐方式
          legend.box.just = "right",  # 图例框架的对齐方式
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),  # 减少图例周围的边距
          legend.key.size = unit(0.5, "cm")  # 设置图例键的大小
        )
      
    }, res = 96)
    
    output$filtered_hist <- renderPlot({
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
          x = paste0("Bins(binsize=", input$bin_size, ")"),  # X轴名称
          y = "Gene Counts",          # Y轴名称
          title = "Histogram of Gene Expression"
        ) +
        guides(color = guide_legend(nrow = 6, byrow = TRUE, title.position = "top")) +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title.x = element_text(vjust = 0.5, hjust = 0.5),  # X轴标题垂直居中
          axis.title.y = element_text(hjust = 0.5),  # Y轴标题水平居中
          panel.border = element_rect(colour = "black", fill = NA, size = 1),  # 给图表加上边框
          legend.box.background = element_rect(colour = "black", fill = NA, size = 0.2),  # 给图例加上边框
          legend.title = element_blank(),
          legend.position = c(1, 1),  # 将图例放在图表内部，位置根据需要调整
          legend.justification = c("right", "top"),  # 图例的对齐方式
          legend.box.just = "right",  # 图例框架的对齐方式
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),  # 减少图例周围的边距
          legend.key.size = unit(0.5, "cm")  # 设置图例键的大小
        )
      
    }, res = 96)
  
    
#    observeEvent(data(), {
#      x <- log2(data()$`b05 Neutrophil1 Exp2` + 1)
#      updateSliderInput(inputId = "bin_size", min = min(x), max = max(x) %/% 2)
#    })  
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
        sliderInput(NS(id, "bin_size"), "Bin size", min = 0, max = 1, value = 0.25, step = 0.25),
        plotOutput(NS(id, "hist"))
      ),
      column(6, 
        plotOutput(NS(id, "grouped_hist"))
      )
    )
  )
}