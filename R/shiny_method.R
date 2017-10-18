plot.interactive <- function(object){
  k.max <- object$tree.max
  k.min <- object$tree.min
  
  non.overlap <- object$non.overlap.clu
  data.tsne <- object$tsne.y
  
  overlap.clu <- object$overlap.clu
  data.pre <- object$data.pre
  
  shinyApp(
    ui = fluidPage(
      
      pageWithSidebar(
        
        #  Application title
        headerPanel("Boosa overlap clustering result"),
        
        # Sidebar with sliders that demonstrate various available options
        sidebarPanel(
          
          sliderInput("hc.k", "K in HC (with recommended choice of K:)", 
                      min = k.min, max = k.max, value = k.min, step = 1),
          sliderInput("overlap.k", "K in Overlap (with recommended choice of K:)", 
                      min = k.min, max = k.max, value = k.min, step = 1),
          selectInput("hc.method", "HC.Method:",
                      list("ward.D" = "ward.D", 
                           "ward.D2" = "ward.D2", 
                           "single" = "single", 
                           "complete" = "complete", 
                           "average" = "average"))
          
        ),
        
        # Show a table summarizing the values entered
        mainPanel(
          #h4("Summary"),
          #verbatimTextOutput("summary")
          tabsetPanel(
            tabPanel("Hierarchical Clustering", rbokehOutput("hc.clust", width = 500, height = 500)),
            tabPanel("Overlap Clustering", rbokehOutput("overlap.clust", width = 800, height = 3000))
            #tabPanel("hmh.group", rbokehOutput("group.rbokeh", width = 600, height = 1200))
          )
        )
      )
    ),
    
    server = function(input, output, session) {
      
      overlap.melt <- reactive({
        overlap.k.index <- input$overlap.k - k.min + 1
        overlap.res <- overlap.clu[[overlap.k.index]]
        cell.overlap.size <- apply(overlap.res[,-1], 1, sum)
        overlap.res <- transform(overlap.res, cell.overlap.size = cell.overlap.size,
                                 X1 = data.tsne[,1], X2 = data.tsne[,2])
        overlap.res.melt <- melt(overlap.res, id.vars = c("cell", "cell.overlap.size", "X1", "X2"))
        transform(overlap.res.melt, value = factor(value))
        
      })
      
      hc.res <- reactive({
        hc.k.index <- input$hc.k - k.min + 1
        as.factor(non.overlap[,hc.k.index])
      })
      
      output$hc.clust <- renderRbokeh({
        
        figure(title = "visualization with tsne", xlab = "tsne.x", ylab = "tsne.y", 
               legend_location = NULL) %>%
          ly_points(x = X1, y = X2, color = hc.res(),
                    data = data.tsne, hover = list(cell)) 
      })
      
      output$overlap.clust <- renderRbokeh({
        n.melt <- dim(overlap.melt())[1]
        idx <- split(1:n.melt, overlap.melt()$variable)
        figs <- lapply(idx, function(x) {
          figure(width = 400, height = 400, legend_location = NULL) %>%
            ly_points(x = X1, y = X2, color = value,
                      data = overlap.melt()[x,], hover = list(cell = paste(":", cell, sep = ""), 
                                                              overlap.size = paste(":", cell.overlap.size, sep = "")), 
                      size = cell.overlap.size*6) 
        })
        
        grid_plot(figs, ncol = 2, same_axes = TRUE, link_data = TRUE)
        
      })
      
    }
  )
}
