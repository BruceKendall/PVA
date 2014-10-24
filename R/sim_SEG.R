sim_SEG <- function() {
  require(shiny)
  require(PVA)
  shinyApp(
    ui = fluidPage(
      # Application title
      titlePanel("Stochastic exponential growth simulation"),
      
      # Sidebar with a slider input for the number of bins
      sidebarLayout(
        sidebarPanel(
          numericInput("mu",
                       "Mean log growth rate (mu)",
                       value=0.01,
                       min=-1,
                       max=1),
          numericInput("sigma2",
                       "Variance of log growth rate (sigma2)",
                       value=0.01,
                       min=0,
                       max=1),
          numericInput("N0",
                       "Initial population size (N0)",
                       value=100,
                       min=1),
          numericInput("time_horizon",
                       "Number of years to simulate",
                       value=10,
                       min=2),
          numericInput("num_sims",
                       "Number of simulations",
                       value=1,
                       min=1),
          br(),
          actionButton("runAgain",
                       "Run Again")
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
          fluidRow(
            column(8,
                   h3("Time series"),
                   plotOutput("tsPlot")
                   
            ),
            column(4,
                   h3("Final distribution"),
                   plotOutput("finalHist"))
          ),
          fluidRow(
            column(8,
                   plotOutput("tsPlotlog")
                   
            ),
            column(4,
                   plotOutput("finalHistlog"))
          )
        )
      )
    ),
    server = function(input,output) {
      Years <- reactive({ 0:input$time_horizon })
      Nt1 <- reactive({
        input$runAgain
        replicate(input$num_sims, 
                  c(input$N0, 
                    simulateSEG(input$mu, input$sigma2, input$N0, input$time_horizon))
        )
      })
      
      
      
      output$tsPlot <- renderPlot({
        matplot(Years(), Nt1(), type='l', xlab="Time (t)", ylab=expression(N[t])) 
      })
      output$tsPlotlog <- renderPlot({
        matplot(Years(), Nt1(), type="l", log="y", xlab="Time (t)", 
                ylab=expression(N[t]~~(log~~scale)))
      })
      output$finalHist <- renderPlot({
        hist(as.matrix(Nt1()[(input$time_horizon+1),]),
             main="", xlab=expression(N[t]))
      })
      output$finalHistlog <- renderPlot({
        hist(as.matrix(log(Nt1()[(input$time_horizon+1),])),
             main="", xlab=expression(log(N[t])))        
      })
      
      
    }
  )
}