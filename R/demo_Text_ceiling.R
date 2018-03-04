demo_Text_ceiling <- function() {
  require(shiny)
  require(PVA)
  shinyApp(
    ui = fluidPage(
      # Application title
      titlePanel("Mean time to extinction with ceiling density dependence"),
      
      # Sidebar with a slider input for the number of bins
      sidebarLayout(
        sidebarPanel(
          sliderInput("mu",
                       "Low-density growth rate (mu)",
                       value=0.01,
                      min=-0.5, max=0.5,
                       step=0.01),
          sliderInput("sigma2",
                       "Variance of log growth rate (sigma2)",
                       value=0.05,
                       min=0,
                       max=1,
                       step=0.01),
          sliderInput("Kmax",
                      "Maximum K on plot",
                      value=100,
                      min=10,
                      max=1000,
                      step=10)
        ),
        
        # Show the plot
        mainPanel(
          plotOutput("TbarPlot")
        )
      )
    ),
    server = function(input,output) {
      Tbar <- reactive({
        Kvec <- 2:input$Kmax
        c <- input$mu/input$sigma2
        (Kvec^(2*c) - 1 - 2*c*log(Kvec))/(2*input$mu*c)
      })
      # Plot the result
      output$TbarPlot <- renderPlot({
        plot(2:input$Kmax, Tbar(), type='l', ylim=c(0,500), xlab="K", 
             ylab="Mean time to extinction") 
     })
      
      
    }
  )
}