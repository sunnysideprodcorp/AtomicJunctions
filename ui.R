shinyUI(fluidPage(theme="style.css",

    headerPanel("Atomic point contact break-junction traces"),

    singleton(tags$head(HTML(
      '
    <script type="text/javascript">
    $(document).ready(function() {
      
      // disable start_proc button after a click
      Shiny.addCustomMessageHandler("disableButton", function(message) {
      $("#clusterButton").attr("disabled", "true");
      });
      
      // Enable start_proc button when computation is finished
      Shiny.addCustomMessageHandler("enableButton", function(message) {
      $("#clusterButton").removeAttr("disabled");
      });

      // disable start_proc button after a click
      Shiny.addCustomMessageHandler("disableButton2", function(message) {
      $("#button2D").attr("disabled", "true");
      });
      
      // Enable start_proc button when computation is finished
      Shiny.addCustomMessageHandler("enableButton2", function(message) {
      $("#button2D").removeAttr("disabled");
      });

      })
      </script>
    '
    ))),

    #row 1
    fluidRow(
      column(1),
      column(3, align = "center", wellPanel(
        h4("See a sample trace"),
          fluidRow( column(6, textInput("startGsteps", label = "high G", value = "1.5")),
                    column(6, textInput("stopGsteps", label = "low G", value = ".3"))
            
          ),
        sliderInput("sensitivitySlider", label = "Adjust sensitivity",
                    min = 20, max = 200, value = 100),
        actionButton('samp','New Sample' )
      )),
      column(3,
        plotOutput("trace", height="300px", width = "300px")
      ),
      column(3,
        plotOutput("stepsShow", height="300px", width = "300px")
      )
      
    ),
    
    #row 2 - 2d histogram
    fluidRow(
      column(1),
      column(3, align = "center", wellPanel(
        h4("Two dimensional histogram"),
        fluidRow(column(4), column(4, textInput("align2D", label = "align at", value = ".3"))),
        actionButton('button2D', '2D Histogram')

      )),
      column(2),
      column(3,
             plotOutput("hist2D", height="300px", width = "350px")
      )

      ),
    #row 3 - clustering
    fluidRow(
      column(1), 
      column(3, align = "center", wellPanel(        h4("Within-plateau step clustering"),
                                                    p("Click the button to see clusters of steps identified"),
                                                    sliderInput("clusterSlider", label = "Set number of clusters",
                                                                min = 2, max = 8, value = 3),
                                                    actionButton('clusterButton', 'See Clusters')
      )),
      column(2),
      column(3, plotOutput("clustersShow", height="300px", width = "350px"))
    )
    

    

    
    
    ))
