cdf_SEG <- function() {
  require(shiny)
  require(PVA)
  shinyApp(
    ui = fluidPage(
      # Application title
      titlePanel("Cumulative quasi-extinction risk under stochastic exponential growth"),
      
      # Sidebar with a slider input for the number of bins
      sidebarLayout(
        sidebarPanel(
          numericInput("mu",
                       "Mean log growth rate (mu)",
                       value=-0.05,
                       min=-1,
                       max=1,
                       step=0.01),
          numericInput("sigma2",
                       "Variance of log growth rate (sigma2)",
                       value=0.03,
                       min=0,
                       max=1,
                       step=0.01),
          numericInput("N0",
                       "Initial population size (N0)",
                       value=100,
                       min=1),
          numericInput("Nx",
                       "Quasi-extinction threshold (Nx)",
                       value=50,
                       min=1),
          numericInput("time_horizon",
                       "Number of years to simulate",
                       value=100,
                       min=2),
          numericInput("num_sims",
                       "Number of simulations",
                       value=1000,
                       min=2)
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
          plotOutput("cdfPlot")
#          textOutput("testText")
        )
      )
    ),
    server = function(input,output) {
      Years <- reactive({ 0:input$time_horizon })
      d <- reactive({ log(input$N0/input$Nx) })
      epsilons <- reactive({ 
        nn <- input$num_sims * input$time_horizon
        matrix(data = rnorm(n = nn), nrow = input$time_horizon)
      })
      mus <- reactive({ 
#        seq(from=input$mu, by=input$mu, length=input$time_horizon) 
        input$mu * (1:input$time_horizon)
      })
      mins <- reactive({ 
        apply(epsilons()*sqrt(input$sigma2) + mus(), 2, cummin)
      })
      CDF <- reactive({
        exts <- mins() < -d()
        c(0, apply(exts, 1, mean))
      })
      
#      output$testText <- renderText({d()})
      output$cdfPlot <- renderPlot({
        plot(Years(), CDF(), type='l', xlab="Years", 
             ylab="Cumulative extinction risk") 
     })
      
      
    }
  )
}