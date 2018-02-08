#' Extinction risk demo: Ricker model
#' 
#' A Shiny demo of extinction risk under the SEG model. The user can set the low-density growth rate,
#' the equilibrium abundance,
#' the variance of the log growth rate, the initial population size, the quasi-extinction population size,
#' and the length of simulation.
#' The number of replicate simulations can also be changed.
#'
#' A single graph is produced, showing the cumulative disribution function (CDF) 
#' of (quasi-) extinction risk. This is calculated using \code{\link{extRiskRicker}}.
#' 
#' The simulations update immediately when you change any parameter value. 
#' 
#' 
#' @return Nothing is returned; the commmand is run for its side effect of 
#' launching a shiny app. You will have to close the shiny window or hit the
#' Stop button in the console to get the console prompt back.
#' 
demo_cdf_Ricker <- function() {
  require(shiny)
  require(PVA)
  shinyApp(
    ui = fluidPage(
      # Application title
      titlePanel("Cumulative quasi-extinction risk under Ricker density dependence"),
      
      # Sidebar with a slider input for the number of bins
      sidebarLayout(
        sidebarPanel(
          numericInput("r",
                       "Low-density growth rate (r)",
                       value=1.5,
                       min=0,
                       step=0.1),
          numericInput("K",
                       "Equilibrium abundance (K)",
                       value=100,
                       min=0),
          numericInput("sigma2",
                       "Variance of log growth rate (sigma2)",
                       value=0.05,
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
      CDF <- reactive({
        extRiskRicker(input$r, input$K, input$sigma2, input$time_horizon, input$N0, 
                      input$num_sims, input$Nx) 
      })
      # Plot the CDF
      output$cdfPlot <- renderPlot({
        plot(CDF~Time, data=CDF(), type='l', xlab="Years", 
             ylab="Cumulative extinction risk", ylim = c(0, 1)) 
     })
      
      
    }
  )
}