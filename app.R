library(shiny)

# 加载UI定义
source('ui.R')
# 加载服务器逻辑
source('server.R')

# 启动应用
shinyApp(ui = ui, server = server)