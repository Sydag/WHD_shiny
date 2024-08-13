library(shiny)

# logifySlider javascript function
JS.logify <-
  "
// function to logify a sliderInput
function logifySlider (sliderId, sci = false) {
  if (sci) {
    // scientific style
    $('#'+sliderId).data('ionRangeSlider').update({
      'prettify': function (num) { return ('10<sup>'+num+'</sup>'); }
    })
  } else {
    // regular number style
    $('#'+sliderId).data('ionRangeSlider').update({
      'prettify': function (num) { return (Math.pow(10, num)); }
    })
  }
}"

# call logifySlider for each relevant sliderInput
JS.onload <-
  "
// execute upon document loading
$(document).ready(function() {
  // wait a few ms to allow other scripts to execute
  setTimeout(function() {
    // include call for each slider
    logifySlider('log_slider', sci = false)
    logifySlider('log_slider2', sci = true)
  }, 5)})
"

ui <- fluidPage(
  tags$head(tags$script(HTML(JS.logify))),
  tags$head(tags$script(HTML(JS.onload))),
  
  sliderInput("log_slider", "Log Slider (numbers):",
              min = -5, max = 3, value = -4, step = 1),
  
  sliderInput("log_slider2", "Log Slider (sci. notation):",
              min = -8, max = -1, value = -6, step = 1),
  
  br(),
  
  textOutput("readout1"),
  textOutput("readout2")
)

server <- function(input, output, session) {
  output$readout1 <- reactive({
    paste0("Selected value (numbers): ", input$log_slider, " = ", 10^input$log_slider)
  })
  
  output$readout2 <- reactive({
    paste0("Selected value (sci. notation): ", input$log_slider2, " = ", 10^input$log_slider2)
  })
}

shinyApp(ui, server)
