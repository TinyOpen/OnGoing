
#' Opens \code{BOSSA} results in an interactive session in a web browser.
#'
#' Runs interactive \code{shiny} session of \code{BOSSA} based on precomputed clusterings.
#'
#' @param object a output of \code{BossaClust} function
#'
#' @return Opens a browser window with an interactive \code{shiny} app and visualize
#' all precomputed clusterings.
#'
#' @name bossa_interactive
#'
#' @importFrom shiny fluidPage pageWithSidebar mainPanel headerPanel pageWithSidebar selectInput shinyApp sliderInput tabPanel tabsetPanel  pageWithSidebar reactive
#' @importFrom rbokeh renderRbokeh ly_points figure grid_plot
#' @importFrom reshape2 melt
#' @export
#'

bossa_interactive <- function(object){
  k.max <- object$tree.max
  k.min <- object$tree.min
  
  non.overlap <- object$non.overlap.clu
  data.tsne <- as.data.frame(object$tsne.y)
  colnames(data.tsne) <- c("X1", "X2", "cell")
  
  overlap.clu <- object$overlap.clu
  non.overlap.clu <- object$non.overlap.clu
  
  de.p.overlap <- object$de.p.overlap
  de.p.non.overlap <- object$de.p.non.overlap 
  
  
  hc.de.plot <- object$hc.de.plot
  overlap.de.plot <- object$overlap.de.plot
  bef.de.plot <- object$bef.de.plot
  
  overlap.melt.data.ori <- object$overlap.melt.data
  
  max.sub.k <- dim(overlap.clu[[1]])[2] - 1
  
  ori.overlap <- as.data.frame(object$ori.overlap[,-c(1,2)])
  mer.clu <- object$mer.clu
  if(length(mer.clu) == 0) {
    mer.clu.zero <- TRUE
    mer.clu <- list(clus.warning = "there is no cluster to be merged")
  } else {
    mer.clu.zero <- FALSE
    mer.clu <- as.list(object$mer.clu)
  }
  
  

  shinyApp(
    
    ui = navbarPage(
      "BoosaClust",
      tabPanel("HC & heatmap",
               fluidPage(
                 fluidRow(
                   column(2,
                          selectInput(inputId = "hc.k",
                                      label = "K in HC:",
                                      choices = as.list(seq(k.min, k.max, 1)),
                                      selected = "k.min")
                   ),
                   column(10,
                          tabsetPanel(
                            tabPanel(
                              "HC tsne plot",
                              plotOutput("hc.clust")
                            ),
                            tabPanel(
                              "DE genes",
                              fluidRow(
                                column(8,
                                       h4(textOutput("hint"))
                                       )
                              ),
                              fluidRow(
                                column(10,
                                       d3heatmapOutput("heatmap.hc", height = 600)
                                ),
                                column(2,
                                       DT::dataTableOutput("hc.de.table")
                                )
                              )
                            )
                          )
                   )
                 )
                 )
               ), 
                 
      
      tabPanel("Overlap Clust & heatmap", 
               sidebarPanel(width = 2,
                            selectInput(inputId = "overlap.k",
                                        label = "K in Overlap:",
                                        choices = as.list(seq(k.min, k.max, 1)),
                                        selected = "k.min")
               ),
               mainPanel(width = 10,
                         tabsetPanel(
                           tabPanel("Overlap tsne plot",
                                    plotOutput("overlap.clust")),
                           tabPanel("DE genes",
                             fluidPage(
                               fluidRow(
                                 column(10,
                                        selectInput("overlap.sub.k", 
                                                    label = paste("There are", max.sub.k, "subclusters"),
                                                    choices = as.list(seq(1, max.sub.k, 1)),
                                                    selected = "1"),
                                        d3heatmapOutput("heatmap.overlap", height = 600)
                                        
                                 ),
                               fluidRow(
                                 column(2,
                                        DT::dataTableOutput("over.de.table"))
                               )
                             )
                           )
                         )
                         )
               )
      ),
      if(!mer.clu.zero){
      tabPanel("DE Before Merge", 
               fluidPage(
                 fluidRow(
                   column(7,
                          selectInput(inputId = "mer.clu.name",
                                      label = "find DE before merging:",
                                      choices = mer.clu,
                                      selected = mer.clu[[1]]),
                          d3heatmapOutput("heatmap.bef.merge", height = 600)
                   ),
                   column(2,
                          DT::dataTableOutput("bef.de.table")),
                 column(3,
                        plotOutput("de.bef.clust")
                 )
               )
               )
      )}
      ),
               
              
    server = function(input, output, session) {
      
      
      observe({
        x <- as.numeric(input$overlap.k) - k.min + 1
        y <- dim(overlap.clu[[x]])[2] - 1 
        updateSelectInput(session, inputId = "overlap.sub.k", 
                          label = paste("There are", y, "subclusters.
                                        Choose one to see the heatmap."),
                          choices = as.list(1:y),
                          selected = "1")
        
      })
      
      lap.k <- reactive({
        as.numeric(input$overlap.k) - k.min + 1
      })
      
      hc.k <- reactive({
        kc.new.k <- as.numeric(input$hc.k) - k.min +1
        kc.new.k
      })
      
      overlap.melt <- reactive({
        K = lap.k()
        overlap.melt.data.ori[[K]]
      })
      
      hc.res <- reactive({
        hc.k.index <- hc.k()
        as.factor(non.overlap[,hc.k.index])
      })
      
      mer.ind <- reactive({
        mer.name <- as.character(input$mer.clu.name)
        unlist(strsplit(mer.name, ", ", fixed = TRUE))
      })
      
      heatmap.bef <- reactive({
        mer.de.ind <- which(mer.clu == as.character(input$mer.clu.name))
        bef.de.data <- as.data.frame(t(bef.de.plot[[mer.de.ind]]))
        
      })
      
      # tsne visualization for HC clustering
      output$hc.clust <- renderPlot({
        cluster <- hc.res()
        ggplot(data = data.tsne, aes(x = X1, y = X2))+
          geom_point(aes(color = cluster),
                     alpha = 0.8, size = 3.5 )+
          ggtitle("HC Clustering with Recommended K") +
          labs(x = "tsne.x", y = "tsne.y")+
          theme_minimal() +
          theme(
            legend.title = element_text(size = 12, color = "salmon", face = "bold"),
            legend.justification = c(0, 1),
            legend.position = "right",
            legend.background = element_blank(),
            legend.key = element_blank()
          ) +
          labs(subtitle = "Visualization with tsne")
        
      })
      
      # tsne visualization for Overlap clustering
      output$overlap.clust <- renderPlot({
        overlap.melt.data <- overlap.melt()
        colnames(overlap.melt.data) <- c("cell", "overlap.size", "tsne.x",
                                         "tsne.y", "cluster.level", "belong")
        overlap.melt.data$overlap.size <- as.factor(overlap.melt.data$overlap.size)
        overlap.melt.data.filter <- subset(overlap.melt.data, belong == 1)
        ggplot(data = overlap.melt.data.filter, aes(x = tsne.x, y = tsne.y, group = cluster.level))+
          geom_point(aes(size = overlap.size),
                     alpha = 0.8, size = 3, color = "salmon")+
          
          facet_wrap(~cluster.level, ncol = 2)+
          ggtitle("Overlap Clusters") +
          labs(x = "tsne.x", y = "tsne.y")+
          theme_minimal() +
          theme(
            legend.title = element_text(size = 12, color = "salmon", face = "bold"),
            legend.justification = c(0, 1),
            legend.position = "right",
            legend.background = element_blank(),
            legend.key = element_blank()
          ) +
          labs(subtitle = "Visualization with tsne")
        
      }) 
      
      output$heatmap.hc <- renderD3heatmap({
        hc.de.plot.data <- as.data.frame(t(hc.de.plot[[hc.k()]]))
        my.height <- dim(hc.de.plot.data)[1]
        my.width <- dim(hc.de.plot.data)[2]
        
        d3heatmap(hc.de.plot.data, cexCol = 0.2 + 50/my.height, 
                  cexRow = 0.2 + 50/my.width,
                  quoted = TRUE, dendrogram = "none")
      })
      
      
      heatmap.over <- reactive({
        overlap.de.plot.big <- overlap.de.plot[[lap.k()]]
        overlap.de.plot.small <- overlap.de.plot.big[[as.numeric(input$overlap.sub.k)]]
        as.data.frame(t(overlap.de.plot.small))
      })
      
      output$heatmap.overlap <- renderD3heatmap({
        my.height <- dim(heatmap.over())[1]
        my.width <- dim(heatmap.over())[2]
        d3heatmap(heatmap.over(), dendrogram = "none", 
                  cexCol = 0.2 + 50/my.height, 
                  cexRow = 0.2 + 50/my.width)
      })
      
      
      
      output$de.bef.clust <- renderPlot({
        venn.data <- ori.overlap[, mer.ind()]
        venn(venn.data, ilabels = TRUE, counts = TRUE, zcolor = "style", cexil = 0.8, cexsn = 1)
      })
      
      output$heatmap.bef.merge <- renderD3heatmap({
        d3heatmap(heatmap.bef(), dendrogram = "none")
      })
      
      
      output$hc.de.table <- DT::renderDataTable(DT::datatable({
        hc.de.plot.data <- as.data.frame(t(hc.de.plot[[hc.k()]]))
        gene.name <- rownames(hc.de.plot.data)
        de.len <- length(gene.name)
        hc.de.table <- data.frame(de.gene = gene.name)
        hc.de.table
      }))
      
      output$over.de.table <- DT::renderDataTable(DT::datatable({
        gene.name <- rownames(heatmap.over())
        de.len <- length(gene.name)
        over.de.table <- data.frame(de.gene = gene.name)
        over.de.table
      }))
      
      output$bef.de.table <- DT::renderDataTable(DT::datatable({
        gene.name <- colnames(heatmap.bef())
        de.len <- length(gene.name)
        bef.de.table <- data.frame(de.gene = gene.name)
        bef.de.table
      }))

      output$hint <- renderText({
        hint <- "And you can choose 'K in HC' for another heatmap"
        hint
      })
      
    }
  )
}