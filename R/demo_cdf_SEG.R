#' Extinction risk demo: Stochastic Exponential Growth model
#' 
#' A Shiny demo of extinction risk under the SEG model. The user can set the mean log growth rate,
#' the variance of the log growth rate, the initial population size, the quasi-extinction population size,
#' and the length of simulation.
#' The number of replicate simulations can also be changed.
#'
#' A single graph is produced, showing the cumulative disribution function (CDF) of (quasi-) extinction risk.
#' 
#' The simulations update immediately when you change any parameter value. 
#' The updates are fast because new random numbers are calculated only if the simulation length or number of replicates is changed.
#' 
#' 
#' @return Nothing is returned; the commmand is run for its side effect of 
#' launching a shiny app. You will have to close the shiny window or hit the
#' Stop button in the console to get the console prompt back.
#' 
demo_cdf_SEG <- function() {
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
        
        # Show the plot
        mainPanel(
          plotOutput("cdfPlot")
        )
      )
    ),
    server = function(input,output) {
      # Note that the following code could have been replaced with a simple
      # call to extRiskSEG(), but for some parameters this longer form will
      # update faster when the parameter is changed
      Years <- reactive({ 0:input$time_horizon })
      d <- reactive({ log(input$N0/input$Nx) })
      
      # Random numbers are only updated if number of simulations or time 
      # horizon is changed
      epsilons <- reactive({ 
        nn <- input$num_sims * input$time_horizon
        apply(matrix(data = rnorm(n = nn), nrow = input$time_horizon), 2, cumsum)
      })
      
      # Vector of expected changes from zero
      mus <- reactive({ 
        input$mu * (1:input$time_horizon)
      })
  
      # This both converts to actual log values (relative to zero)
      # and takes the cumulative minimum. It's the only place that
      # depends on sigma2
      mins <- reactive({ 
        apply(epsilons()*sqrt(input$sigma2) + mus(), 2, cummin)
      })
      
      # And, finally, calculate the CDF
      CDF <- reactive({
        exts <- mins() < -d()
        c(0, apply(exts, 1, mean))
      })
      
      # Plot the CDF
      output$cdfPlot <- renderPlot({
        plot(Years(), CDF(), type='l', xlab="Years", 
             ylab="Cumulative extinction risk") 
     })
      
      
    }
  )
}