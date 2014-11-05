sim_Ricker <- function() {
  require(shiny)
  require(PVA)
  shinyApp(
    ui = fluidPage(
      # Application title
      titlePanel("Stochastic Ricker model simulation"),
      
      # Sidebar with a slider input for the number of bins
      sidebarLayout(
        sidebarPanel(
          numericInput("r",
                       "Low-density growth rate (r)",
                       value=1,
                       min=0),
          numericInput("K",
                       "Equilibrium abundance (K)",
                       value=100,
                       min=0),
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
                   plotOutput("tsPlot")
                   
            ),
            column(4,
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
        simulateRicker(input$r, input$K, input$sigma2, input$N0, 
                    input$time_horizon, input$num_sims, TRUE)
        
      })
      
      
      
      output$tsPlot <- renderPlot({
        matplot(Years(), Nt1(), type='l', xlab="Time (t)", ylab=expression(N[t]),
                main="Time series") 
      })
      output$tsPlotlog <- renderPlot({
        matplot(Years(), Nt1(), type="l", log="y", xlab="Time (t)", 
                ylab=expression(N[t]~~(log~~scale)))
      })
      output$finalHist <- renderPlot({
        hist(as.matrix(Nt1()[(input$time_horizon+1),]),
             main="Final Distribution", xlab=expression(N[t]))
      })
      output$finalHistlog <- renderPlot({
        hist(as.matrix(log(Nt1()[(input$time_horizon+1),])),
             main="", xlab=expression(log(N[t])))        
      })
      
      
    }
  )
}