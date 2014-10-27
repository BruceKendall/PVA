cdfMap_SEG <- function() {
  require(shiny)
  require(PVA)
  shinyApp(
    ui = fluidPage(
      # Application title
      titlePanel("Cumulative quasi-extinction risk under stochastic exponential growth"),
      
      # Sidebar with a slider input for the number of bins
      sidebarLayout(
        sidebarPanel(
          sliderInput("muRange",
                       "Range for mean log growth rate (mu)",
                       value=c(-0.1, 0.1),
                       min=-1,
                       max=1,
                       step=0.01),
          sliderInput("sigma2Range",
                       "Range for variance of log growth rate (sigma2)",
                       value=c(0, 0.07),
                       min=0,
                       max=1,
                       step=0.001),
          numericInput("nstep",
                      "Number of values on each axis",
                      value=21,
                      min=5,
                      max=100),
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
                       value=10,
                       min=2),
          numericInput("num_sims",
                       "Number of simulations",
                       value=1000,
                       min=2),
          submitButton("Update now")
        ),
        
        # Show the plot
        mainPanel(
          plotOutput("cdfContourPlot")
        )
      )
    ),
    server = function(input,output) {
      # Note that the following code could have been replaced with a simple
      # call to extRiskSEG(), but for some parameters this longer form will
      # update faster when the parameter is changed
      Years <- reactive({ 1:input$time_horizon })
      d <- reactive({ log(input$N0/input$Nx) })
      muList <- reactive({ 
        seq(input$muRange[1], input$muRange[2], length=input$nstep)
      })
      sigma2List <- reactive({ 
        seq(input$sigma2Range[1], input$sigma2Range[2], length=input$nstep)
      })       
      
      # Random numbers are only updated if number of simulations or time 
      # horizon is changed
      epsilons <- reactive({ 
        nn <- input$num_sims * input$time_horizon
        matrix(data = rnorm(n = nn), nrow = input$time_horizon)
      })
      
      CDFs <- reactive({
        # A container for the results
        zz <- matrix(0, input$nstep, input$nstep)
        
        # Cycle through all the combos
        for (i in 1:input$nstep) { # mu loop
          mu <- muList()[i]
          means <- -d() - mu*Years()
          for (j in 1:input$nstep) { # sigma2 loop
            sigma <- sqrt(sigma2List()[j])
            if (sigma == 0) sigma <- 1e-6
            Ext <- epsilons() < (means/sigma)
            zz[i,j] <- sum(apply(Ext, 2, max))
          }
        }
        
        t(zz)/input$num_sims
      })
     
      
      # Make the contour plot
      output$cdfContourPlot <- renderPlot({
        contour(sigma2List(), muList(), CDFs(),
                xlab=expression(sigma^2), ylab=expression(mu),
                main="quasi-extinction probability",
                levels=c((1:9)/1000, (1:9)/100, (1:10)/10))
        abline(h=0, lty=2)
      })
      
      
    }
  )
}