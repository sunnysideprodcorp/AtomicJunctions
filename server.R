source("makeTraces.R")
source("analysis.R")
source("clusters.R")
library(zoo)
library(KernSmooth)
library(ggplot2)

traceNew = makeTrace()

shinyServer(function(input, output, session) {
  
  # initial loading of plots:
  output$trace <- renderPlot({
    plot(traceNew, type='l')
  })
  output$stepsShow <- renderPlot({
    findJumps(beginG =as.numeric(input$startGsteps), endG = as.numeric(input$stopGsteps), trace = traceNew, sensitivityDivisor = input$sensitivitySlider, loq=1)
  })
  

  #row 1
      observe({
        if(input$samp > 0) {
          traceNew <<- makeTrace()
          output$trace <- renderPlot({
            plot(traceNew, type='l')
          })
          
          output$stepsShow <- renderPlot({
              findJumps(beginG = as.numeric(input$startGsteps), endG = as.numeric(input$stopGsteps), trace = traceNew, sensitivityDivisor = input$sensitivitySlider, loq=1)
          })
        }
      })
      
      observe({
        if(input$sensitivitySlider > 0 || (!is.na(as.numeric(input$startGsteps)) && as.numeric(input$startGsteps) > 0) || (!is.na(as.numeric(input$stopGsteps)) && as.numeric(input$stopGsteps) > 0)){
          output$stepsShow <- renderPlot({
                  findJumps(beginG = as.numeric(input$startGsteps), endG = as.numeric(input$stopGsteps), trace = traceNew, sensitivityDivisor = input$sensitivitySlider, loq=1)
          })
        }
      })



  #row 2 - 2d histogram
    observe({
    if(input$button2D > 0) {
     
      output$hist2D <- renderPlot({
        session$sendCustomMessage("disableButton2","button2D")
        # Create a Progress object
        progress2 <- shiny::Progress$new( min = 0, max = 1)
        progress2$set(message = "Computing data", value = 0)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress2$close())
        
        updateProgress <- function(value = NULL, detail = NULL) {
          if (is.null(value)) {
            value <- progress2$getValue()
          }
          progress2$set(value = value, detail = detail)
        }
        h = make2DHist(alignG = as.numeric(input$align2D), updateProgress = updateProgress)
        session$sendCustomMessage("enableButton2","button2D")
        h
      })
    }
  })
  
  
  #row 3 - clustering plot

     output$clustersShow <- renderPlot({
       if(input$clusterButton > 0) {
         # Create a Progress object
         progress <- shiny::Progress$new( min = 0, max = 1)
         progress$set(message = "Computing data", value = 0)
         # Close the progress when this reactive exits (even if there's an error)
         on.exit(progress$close())
         
         updateProgress <- function(value = NULL, detail = NULL) {
           if (is.null(value)) {
             value <- progress$getValue()
           }
           progress$set(value = value, detail = detail)
         }
         
         session$sendCustomMessage("disableButton","clusterButton")
         clust.size <- input$clusterSlider
         a = eval(clust.size)
         numTraces <-30
         list_df <- makeClusters(numTraces, updateProgress = updateProgress)
         df1     <- list_df[[1]]
         drops     <- c("Trace")
         df1.clust <- df1[,!(names(df1) %in% drops)]
         df2            <- bestCluster(df = df1.clust, k = clust.size, df2 = df1, updateProgress = updateProgress)
         df2$clusterCat <- factor(df2$Cluster)
         session$sendCustomMessage("enableButton","clusterButton")
         sizeToUse = eval(length(levels(factor(df2$Cluster))))
         ggplot(df2, aes(x=Noise, y = AvgC), colour=clusterCat) + geom_point( aes(color=clusterCat)) 

       }
     })
})