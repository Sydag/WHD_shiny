library(shiny)
library(plotly)
library(xtable)



#source("/home/pier/Documents/Research/WHD_model/model/WHD.R")
#setwd("/home/pier/Documents/Research/WHD_model/WHD_app/")
source("WHD.R")
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
    logifySlider('x_ini1', sci = true)
    logifySlider('x_ini2', sci = true)
    logifySlider('mu_x1', sci = true)
    logifySlider('mu_x2', sci = true)
  }, 5)})
"

# Define UI : we're using fluidpage to control the layout.
ui <- navbarPage("WHD model", #inverse=TRUE,
                 theme = bslib::bs_theme(bootswatch = "darkly"),
                 tabPanel("System definition",
                          fluidPage(
                            fluidRow(
                              column(8,
                                     img(width="100%", src="fig_schema_model_reversed_new.png")
                                     #uiOutput("pdf"),
                                     #imageOutput("png")
                              ),
                              column(4,
                                     uiOutput("eq_WHD")
                              ),
                              # column(2,
                              #        h3("Display options"),
                              #        wellPanel(
                              #        checkboxInput("dark","Dark theme mode",value = TRUE)
                              #        )
                              # )
                            )
                          )
                          
                          #uiOutput("pdfview")
                 ),
                 tabPanel("Genotype parameters",
                          fluidPage(
                            #ui pars
                            # Application title
                            #titlePanel("WHD_model"),
                            
                            # Custom made layout (not using sidepanel/mainpanel)
                            fluidRow(
                              column(5,
                                     h4("Genotype 1"),
                                     uiOutput("pars_G1")
                                     
                              ),
                              column(5,
                                     h4("Genotype 2"),
                                     uiOutput("pars_G2")
                                     
                              ),
                              column(2,
                                     h4("Variation in parameters"),
                                     uiOutput("var_pars_panel1")
                              )
                            )
                          )
                          
                 ),
                 tabPanel("Deterministic simulations",
                          tabsetPanel(
                            tabPanel("2D plots",
                                     fluidRow(
                                       column(4,
                                              #conditionalPanel("input.Xplot == true",
                                              plotlyOutput("dyn_plotly_x")
                                              #)
                                       ),
                                       column(4,
                                              #conditionalPanel("input.Yplot == true",
                                              plotlyOutput("dyn_plotly_y")      
                                              #)
                                       ),
                                       column(4,
                                              #conditionalPanel("input.Zplot == true",
                                              plotlyOutput("dyn_plotly_z")      
                                              #)
                                       ),
                                     )
                            ),
                            tabPanel("3D plots",
                                     plotlyOutput("phase_plot_3D",height = "900px")
                            )
                          ),
                          #br(),
                          fluidRow(
                            column(3,
                                   uiOutput("cond_G1")
                            ),
                            column(3,
                                   uiOutput("cond_G2")
                            ),
                            column(3,
                                   h4("Joint options"),
                                   wellPanel(
                                     fluidRow(
                                       column(6,
                                              h5("Graphical Opt."),
                                              checkboxInput("clrb","Alt. colour scheme",value = FALSE),
                                              uiOutput("joint_opt"),
                                              actionButton("res_opt","Reset to default",class="btn-dark")
                                              
                                       ),
                                       column(6,
                                              h5("Simulation Opt."),
                                              sliderInput("tmax","Time points", min=100, max=1000, value=600, step=50, ticks=FALSE)
                                              
                                       )
                                     )
                                   )
                            ),
                            column(3,
                                   h4("Additional information"),
                                   uiOutput("splash_dyn"),
                                   h4("Variation in parameters"),
                                   uiOutput("var_pars_panel2")
                            )
                          )
                          
                          
                 ),
                 tabPanel("Stochastic simulations",
                          fluidRow(
                            column(4,
                                   plotlyOutput("stoch_plotly_G1")
                            ),
                            column(4,
                                   plotlyOutput("stoch_plotly_G2")
                            ),
                            column(4,
                                   plotlyOutput("surv_dead")
                            )
                            
                          ),
                          fluidRow(
                            column(4,
                                   fluidRow(
                                     column(8,
                                            uiOutput("stoch_cond_G1")
                                     ),
                                     column(4,
                                            h4("Options G.1"),
                                            wellPanel(
                                              sliderInput("npts1","Number of Ind. per time point", min=2, max=25, value=10, step=1, ticks=FALSE),
                                              uiOutput("scale_opt_stoch1")
                                            )                                          
                                            
                                            #)
                                     )
                                   )
                            ),
                            column(4,
                                   fluidRow(
                                     column(8,
                                            uiOutput("stoch_cond_G2")
                                     ),
                                     column(4,
                                            
                                            h4("Options G.2"),
                                            wellPanel(
                                              sliderInput("npts2","Number of Ind. per time point", min=2, max=25, value=10, step=1, ticks=FALSE),
                                              uiOutput("scale_opt_stoch2")
                                            )
                                     )
                                   )
                            ),
                            
                            column(4,
                                   fluidRow(
                                     column(6,
                                            h4("Survival plot options"),
                                            uiOutput("surv_opt"),
                                     ),
                                     column(6,
                                            h4("Additional information"),
                                            uiOutput("splash_stoch"),
                                     )
                                   ),
                                   h4("Variation in parameters"),
                                   uiOutput("var_pars_panel3")
                                   
                                   #plotlyOutput("hist_dead")
                            )
                          )
                          
                 ),
                 tags$script(HTML("var header = $('.navbar > .container-fluid');
header.append('<div style=\"float:right\"><a href=\"https://edb.cnrs.fr/\"><img src=\"Logo_UT3.jpg\" alt=\"alt\" style=\"float:right;width:160px;height:50px;padding-right:10px;\"> </a></div>');
    console.log(header)")),
                 tags$script(HTML("var header = $('.navbar > .container-fluid');
header.append('<div style=\"float:right\"><a href=\"https://gulbenkian.pt/ciencia/\"><img src=\"logo_Gulbenkian.png\" alt=\"alt\" style=\"float:right;width:160px;height:60px;padding-right:10px;\"> </a></div>');
    console.log(header)")),
                 tags$script(HTML("var header = $('.navbar > .container-fluid');
header.append('<div style=\"float:right\"><a href=\"https://www.ed.ac.uk/biology/\"><img src=\"logo.png\" alt=\"alt\" style=\"float:right;width:200px;height:50px;padding-top:0px;\"> </a></div>');
    console.log(header)"))
                 # tabPanel("Simulations",
                 #          plotlyOutput("dyn_plot_3D",height = "900px")
                 # )
                 
                 
)
# Define server logic required to draw a histogram
server <- function(input, output) {
  thematic::thematic_shiny()
  ### use mathjax to display eq system
  output$eq_WHD <- renderUI({
    withMathJax(paste0('$$
                    \\left \\lbrace 
                    \\begin{array}{rl}
                    \\dfrac{dx}{dt}=&x(1-x) - \\delta xy\\\\[1em]
                    \\dfrac{dy}{dt}=&F(x)G(y,z) - \\varphi y\\\\[1em]
                    \\dfrac{dz}{dt}=&\\omega x + \\eta y - \\xi z& 
                    \\end{array}
                    \\right .
                    $$','with','$$F(x)= \\gamma - \\alpha\\dfrac{x^{u}}{1+\\beta^{v}}$$','and','$$G(y,z)= \\dfrac{1}{1+y^{k}+\\psi z^l}$$')
                
    )
  })
  #observeEvent(input$gene)
  output$pdf <- renderUI({
    tags$iframe(style="height:800px; width:100%", src="fig_schema_model_reversed.pdf")
  })
  output$png <- renderImage({
    img(height="800px", width="800px", src="fig_schema_model_reversed.png")
  },deleteFile = FALSE)
  ###create reactive outputs
  
  color_panel_G1 <- reactive({
    if(input$clrb){
      paste0("background: #2F4562")
    }else{
      paste0("background: #1D7874")
    }
  })
  
  color_panel_G2 <- reactive({
    if(input$clrb){
      paste0("background: #C14C32")
    }else{
      paste0("background: #C39059")
    }
  })
  #hue<-c("#6C5B7B","#F67280")
  #colvec<-c("#F8B195","F67280","#C06C84","#6C5B7B","#355C7D")
  
  output$pars_G1 <- renderUI({
    
    wellPanel(style = color_panel_G1(),
              fluidRow(
                column(6,
                       uiOutput("pars_G1_res")
                ),
                column(6,
                       uiOutput("pars_G1_tol")
                )
              )
    )
  })
  
  output$pars_G2 <- renderUI({
    wellPanel(style = color_panel_G2(),
              fluidRow(
                column(6,
                       uiOutput("pars_G2_res")
                ),
                column(6,
                       uiOutput("pars_G2_tol")
                )
              )
    )
  })
  
  output$pars_G1_res <- renderUI({
    validate(
      need(input$gamma_G1_chckbox != "", "Loading...") # display custom message in need
    )
    
    times <- input$res_G1
    div(id=letters[(times %% length(letters)) + 1],
        withMathJax(),
        # wellPanel(style = color_panel_G1(),
        #           h3("Genotype 1"),
        #           fluidRow(
        #             column(6,
        h5("Resistance parameters"),
        wellPanel(
          sliderInput("delta_1","Pathogen destruction rate (\\(\\delta\\))", min=0.001, max=0.1, value=0.05, step=0.001, ticks=FALSE),
          sliderInput("phi_1","Defense decay rate (\\(\\varphi\\))", min=0.001, max=0.1, value=0.01, step=0.001, ticks=FALSE),
          conditionalPanel(
            condition = "input.gamma_G1_chckbox == false",
            sliderInput("gamma_1","Constitutive immunity (\\(\\gamma\\))", min=0, max=20,
                        if(input$gamma_G1_chckbox==FALSE){value=1.8}else{value=0}, 
                        step=0.1, ticks=FALSE),    
          ),
          sliderInput("alpha_1","Immune activation (\\(\\alpha\\))", min=0, max=50, value=12, ticks=FALSE),
          sliderInput("beta_1","Immune saturation (\\(\\beta\\))", min=0, max=5, value=1, step=0.1, ticks=FALSE),
          sliderInput("u_1","Activation shape parameter (\\(u\\))", min=0, max=2, value=0.1, step=0.1, ticks=FALSE),
          sliderInput("v_1","Saturation shape parameter (\\(v\\))", min=0, max=2, value=0.1, step=0.1, ticks=FALSE),
          sliderInput("w_1","Down-regulation shape parameter (\\(k\\))", min=0, max=2, value=0.5, step=0.1, ticks=FALSE),
          conditionalPanel(
            condition = "input.psi_G1_chckbox == false",
            sliderInput("psi_1", "Strength of damage effect (\\(\\psi\\))", min=0, max=10, 
                        if(input$psi_G1_chckbox==FALSE){value=2}else{value=0},
                        step=.1, ticks=FALSE)
          ),
          sliderInput("l_1","Damage effect shape parameter (\\(l\\))", min=0, max=2, value=0.75, step=0.01, ticks=FALSE)
        )
    )
  })
  
  output$pars_G1_tol <- renderUI({
    times <- input$res_G1
    div(id=letters[(times %% length(letters)) + 1],
        withMathJax(),
        #column(6,
        h5("Tolerance parameters"),
        wellPanel(
          sliderInput("omega_1","Damage production rate from pathogens (\\(\\omega\\))", min=0, max=50, value=10, ticks=FALSE),
          sliderInput("eta_1","Damage production rate from defense (\\(\\eta\\))", min=0, max=2, value=0.01, step=0.01, ticks=FALSE),
          sliderInput("ksi_1","Damage repair rate (\\(\\xi\\))", min=0.01, max=2, value=0.01, step=0.01, ticks=FALSE),
          sliderInput("zd_1","Damage threshold (zd)", min=20, max=1000, value=500, step=10, ticks=FALSE)
        ),
        h5("Model functions"),
        wellPanel(
          checkboxInput("psi_G1_chckbox","Disable damage effect on defense prod. (G.1)",value=FALSE),
          checkboxInput("gamma_G1_chckbox","Disable constitutive immunity (G.1)",value=FALSE)
        ),
        br(),
        actionButton("res_G1","Reset G.1 input state",class=ifelse(input$clrb,"btn-warning","btn-primary"))
    )
    
  })
  
  output$pars_G2_res <- renderUI({
    validate(
      need(input$gamma_G2_chckbox != "", "Loading...") # display custom message in need
    )
    times <- input$res_G2
    div(id=letters[(times %% length(letters)) + 1],
        withMathJax(),
        h5("Resistance parameters"),
        wellPanel(
          sliderInput("delta_2","Pathogen destruction rate (\\(\\delta\\))", min=0.001, max=0.1, value=0.05, step=0.001, ticks=FALSE),
          sliderInput("phi_2","Defense decay rate (\\(\\varphi\\))", min=0.001, max=0.1, value=0.01, step=0.001, ticks=FALSE),
          conditionalPanel(
            condition = "input.gamma_G2_chckbox == false",
            sliderInput("gamma_2","Constitutive immunity (\\(\\gamma\\))", min=0, max=20,
                        if(input$gamma_G2_chckbox==FALSE){value=0.5}else{value=0}, 
                        step=0.1, ticks=FALSE),    
          ),
          sliderInput("alpha_2","Immune activation (\\(\\alpha\\))", min=0, max=50, value=12, ticks=FALSE),
          sliderInput("beta_2","Immune saturation (\\(\\beta\\))", min=0, max=5, value=1, step=0.1, ticks=FALSE),
          sliderInput("u_2","Activation shape parameter (\\(u\\))", min=0, max=2, value=0.1, step=0.1, ticks=FALSE),
          sliderInput("v_2","Saturation shape parameter (\\(v\\))", min=0, max=2, value=0.1, step=0.1, ticks=FALSE),
          sliderInput("w_2","Down-regulation shape parameter (\\(k\\))", min=0, max=2, value=0.5, step=0.1, ticks=FALSE),
          conditionalPanel(
            condition = "input.psi_G2_chckbox == false",
            sliderInput("psi_2","Strength of damage effect (\\(\\psi\\))", min=0, max=10,
                        if(input$psi_G2_chckbox==FALSE){value=2}else{value=0},
                        step=.1, ticks=FALSE)
          ),
          sliderInput("l_2","Damage effect shape parameter (\\(l\\))", min=0, max=2, value=0.75, step=0.01, ticks=FALSE)
        )
    )
  })
  
  output$pars_G2_tol <- renderUI({
    times <- input$res_G2
    div(id=letters[(times %% length(letters)) + 1],
        withMathJax(),
        h5("Tolerance parameters"),
        wellPanel(
          sliderInput("omega_2","Damage production rate from pathogens (\\(\\omega\\))", min=0, max=50, value=10, ticks=FALSE),
          sliderInput("eta_2","Damage production rate from defense (\\(\\eta\\))", min=0, max=2, value=0.01, step=0.01, ticks=FALSE),
          sliderInput("ksi_2","Damage repair rate (\\(\\xi\\))", min=0.01, max=2, value=0.01, step=0.01, ticks=FALSE),
          sliderInput("zd_2","Damage threshold (zd)", min=20, max=1000, value=500, step=10, ticks=FALSE)
        ),
        h5("Model functions"),
        wellPanel(
          checkboxInput("psi_G2_chckbox","Disable damage effect on defense prod. (G.2)",value=FALSE),
          checkboxInput("gamma_G2_chckbox","Disable constitutive immunity (G.2)",value=FALSE)
        ),
        br(),
        actionButton("res_G2","Reset G.2 input state",class=ifelse(input$clrb,"btn-warning","btn-primary"))
    )
  })
  
  output$joint_opt <- renderUI({
    times <- input$res_opt
    div(id=letters[(times %% length(letters)) + 1],
        withMathJax(),
        selectInput("scale_dynx","Transform \\(x\\) axis",choices=c("no","log10"),selected = "no"),
        selectInput("scale_dyny","Transform \\(y\\) axis",choices=c("no","log10","log2"),selected = "no")
    )
  })
  
  output$scale_opt_stoch1 <- renderUI({
    times <- input$res_cond1
    div(id=letters[(times %% length(letters)) + 1],
        withMathJax(),
        selectInput("scale_stoch1x","Transform \\(x\\) axis",choices = c("no","log10"),selected = "no"),
        selectInput("scale_stoch1y","Transform \\(y\\) axis",choices = c("no","log10","log2"),selected = "log2"),
        br(),
        actionButton("res_par_stoch1","Reset input state G.1",class=ifelse(input$clrb,"btn-warning","btn-primary"))
    )
  })
  
  output$scale_opt_stoch2 <- renderUI({
    times <- input$res_cond2
    div(id=letters[(times %% length(letters)) + 1],
        withMathJax(),
        selectInput("scale_stoch2x","Transform \\(x\\) axis",choices = c("no","log10"),selected = "no"),
        selectInput("scale_stoch2y","Transform \\(y\\) axis)",choices = c("no","log10","log2"),selected = "log2"),
        br(),
        actionButton("res_par_stoch2","Reset input state G.2",class=ifelse(input$clrb,"btn-warning","btn-primary"))
    )
  })
  
  output$surv_opt <- renderUI({
    times <- input$res_surv
    div(id=letters[(times %% length(letters)) + 1],
        withMathJax(),
        wellPanel(
          sliderInput("nb_surv","Number of simulated host",min=100,max=1000,value=250,step = 10,ticks=FALSE),
          br(),
          actionButton("res_surv","Reset to default",class="btn-dark")
        )
        # selectInput("scale_survx","Transform \\(x\\) axis",choices = c("no","log10"),selected = "no"),
        # selectInput("scale_survy","Transform \\(y\\) axis)",choices = c("no","log10"),selected = "no"),
        #
    )
  })
  
  output$cond_G1 <- renderUI({
    times <- input$res_cond1
    div(id=letters[(times %% length(letters)) + 1],
        tags$head(tags$script(HTML(JS.logify))),
        tags$head(tags$script(HTML(JS.onload))),
        withMathJax(),
        h4("Genotype 1"),
        wellPanel(style = color_panel_G1(),
                  h5("Initial conditions"),
                  wellPanel(
                    sliderInput("x_ini1","Initial pathogen load (\\(x_{0}\\))",min=-8,max=-1,value=-6,ticks=FALSE),
                    numericInput("y_ini1","Immune handicap (reduction in initial defense (\\(\\Delta(y_{0})\\)))",min=0,max=50,value=0,step=1),
                    numericInput("z_ini1","Wound (additional initial damage (\\(\\Delta(z_{0})\\)))",min=0,max=50,value=0,step=1)
                  ),
                  br(),
                  actionButton("res_cond1","Reset input state G.1",class=ifelse(input$clrb,"btn-warning","btn-primary"))
        )
    )
  })
  output$cond_G2 <- renderUI({
    times <- input$res_cond2
    div(id=letters[(times %% length(letters)) + 1],
        tags$head(tags$script(HTML(JS.logify))),
        tags$head(tags$script(HTML(JS.onload))),
        h4("Genotype 2"),
        withMathJax(),
        wellPanel(style = color_panel_G2(),
                  h5("Initial conditions"),
                  wellPanel(
                    sliderInput("x_ini2","Initial pathogen load (\\(x_{0}\\))",min=-8,max=-1,value=-6,ticks=FALSE),
                    numericInput("y_ini2","Immune handicap (reduction in initial defense (\\(\\Delta(y_{0})\\)))",min=0,max=50,value=0,step=1),
                    numericInput("z_ini2","Wound (additional initial damage (\\(\\Delta(z_{0})\\)))",min=0,max=50,value=0,step=1)
                  ),
                  br(),
                  actionButton("res_cond2","Reset input state G.2",class=ifelse(input$clrb,"btn-warning","btn-primary"))
        )
    )
  })
  
  ### Stochastic simulations initial conditions
  output$stoch_cond_G1 <- renderUI({
    times <- input$res_par_stoch1
    div(id=letters[(times %% length(letters)) + 1],
        tags$head(tags$script(HTML(JS.logify))),
        tags$head(tags$script(HTML(JS.onload))),
        withMathJax(),
        h4("Initial conditions G.1"),
        wellPanel(style = color_panel_G1(),
                  h5("Inoculum (initial pathogen load)"),
                  wellPanel(
                    fluidRow(
                      column(6,
                             sliderInput("mu_x1","Mean (\\(\\mu_{x_{0}}\\))",min=-8,max=-1,value=-6,ticks=FALSE)
                      ),
                      column(6,
                             numericInput("var_x1","SD (\\(\\sigma_{log(x_{0})}\\))",min=0,max=0.5,value=0.2,step=0.01)
                      )
                    )
                  ),
                  #br(),
                  h5("Immune handicap (reduction in initial defense)"),
                  wellPanel(
                    fluidRow(
                      column(6,
                             numericInput("mu_y1","Mean (\\(\\Delta(y_{0})\\))",min=0,max=50,value=0,step=1)
                      ),
                      column(6,
                             numericInput("var_y1","SD (\\(\\sigma_{log(\\Delta(y_{0}))}\\))",min=0,max=0.5,value=0,step=0.01)
                      )
                    )
                  ),
                  #br(),
                  h5("Wound (additional initial damage)"),  
                  wellPanel(
                    fluidRow(
                      column(6,
                             numericInput("mu_z1","Mean (\\(\\Delta(z_{0})\\))",min=0,max=50,value=1,step=1)
                      ),
                      column(6,
                             numericInput("var_z1","SD (\\(\\sigma_{log(\\Delta(z_{0}))}\\))",min=0,max=0.5,value=0.2,step=0.01)
                      )
                    )
                  )
        )
    )
  })
  
  output$stoch_cond_G2 <- renderUI({
    times <- input$res_par_stoch2
    div(id=letters[(times %% length(letters)) + 1],
        #tags$head(tags$script(HTML(JS.logify))),
        #tags$head(tags$script(HTML(JS.onload))),
        withMathJax(),
        #fluidRow(
        h4("Initial conditions G.2"),
        wellPanel(style = color_panel_G2(),
                  
                  h5("Inoculum (initial pathogen load)"),
                  wellPanel(
                    fluidRow(
                      column(6,
                             sliderInput("mu_x2","Mean (\\(\\mu_{x_{0}}\\))",min=-8,max=-1,value=-6,ticks=FALSE)
                      ),
                      column(6,
                             numericInput("var_x2","SD (\\(\\sigma_{log(x_{0})}\\))",min=0,max=0.5,value=0.2,step=0.01)
                      )
                    )
                  ),
                  #br(),
                  h5("Immune handicap (reduction in initial defense)"),
                  wellPanel(
                    fluidRow(
                      column(6,
                             numericInput("mu_y2","Mean (\\(\\Delta(y_{0})\\))",min=0,max=50,value=0,step=1)
                      ),
                      column(6,
                             numericInput("var_y2","SD (\\(\\sigma_{log(\\Delta(y_{0}))}\\))",min=0,max=0.5,value=0,step=0.01)
                      )
                    )
                  ),
                  #br(),
                  h5("Wound (additional initial damage)"),  
                  wellPanel(
                    fluidRow(
                      column(6,
                             numericInput("mu_z2","Mean (\\(\\Delta(z_{0})\\))",min=0,max=50,value=1,step=1)
                      ),
                      column(6,
                             numericInput("var_z2","SD (\\(\\sigma_{log(\\Delta(z_{0}))}\\))",min=0,max=0.5,value=0.2,step=0.01)
                      )
                    )
                  )
        )
    )
  })
  ###generate data
  simus <- reactive({
    par_G1 <- as.list(c(delta=input$delta_1,
                        phi=input$phi_1,
                        gamma=input$gamma_1,
                        alpha=input$alpha_1,
                        beta=input$beta_1,
                        u=input$u_1, 
                        v=input$v_1,
                        w=input$w_1, #this is now called k
                        psi=input$psi_1,
                        l=input$l_1,
                        omega=input$omega_1,
                        eta=input$eta_1,
                        ksi=input$ksi_1,
                        zcrit=input$zd_1 #also called zd
    )
    )
    par_G2 <- as.list(c(delta=input$delta_2,
                        phi=input$phi_2,
                        gamma=input$gamma_2,
                        alpha=input$alpha_2,
                        beta=input$beta_2,
                        u=input$u_2, 
                        v=input$v_2,
                        w=input$w_2, #this is now called k
                        psi=input$psi_2,
                        l=input$l_2,
                        omega=input$omega_2,
                        eta=input$eta_2,
                        ksi=input$ksi_2,
                        zcrit=input$zd_2 #also called zd
    )
    )
    
    pdt=seq(0,input$tmax,by=2)
    
    par_ini1=as.list(c(mu_x=10^input$x_ini1,
                       mu_y=20,
                       mu_z=20,
                       sigma_logx=0.15,
                       sigma_logy=0,
                       sigma_logz=0))
    
    par_ini2=as.list(c(mu_x=10^input$x_ini2,
                       mu_y=20,
                       mu_z=20,
                       sigma_logx=0.15,
                       sigma_logy=0,
                       sigma_logz=0))
    
    #readParameters("Documents/Research/WHD_model/model/new_surv_sim/parameters_AC_new_canon_F1.txt")
    #readParameters("Documents/Research/WHD_model/model/new_surv_sim/startingcondition.txt")
    eq_G1 <- equilibrium(par_G1)
    eq_G2 <- equilibrium(par_G2)
    
    ts_G1 <- tsurv(par_ini1[[1]],ifelse(eq_G1[1,2]-input$y_ini1<0,0,eq_G1[1,2]-input$y_ini1),eq_G1[1,3]+input$z_ini1,par_G1,tmax=5000)
    ts_G2 <- tsurv(par_ini2[[1]],ifelse(eq_G2[1,2]-input$y_ini2<0,0,eq_G2[1,2]-input$y_ini2),eq_G2[1,3]+input$z_ini2,par_G2,tmax=5000)
    
    sol_G1 <- dyninf(par_ini1[[1]],ifelse(eq_G1[1,2]-input$y_ini1<0,0,eq_G1[1,2]-input$y_ini1),eq_G1[1,3]+input$z_ini1,pdt,par_G1)
    sol_G2 <- dyninf(par_ini2[[1]],ifelse(eq_G2[1,2]-input$y_ini2<0,0,eq_G2[1,2]-input$y_ini2),eq_G2[1,3]+input$z_ini2,pdt,par_G2)
    
    ###setting up the colour object
    if(input$clrb){
      hue<-c("#2F4562","#C14C32")#,"#F7441B"
    }else{
      hue <- c("#1D7874","#C39059")
      
    }
    
    obj_lst <- return(
      list(
        par_1=par_G1,
        par_2=par_G2,
        eq1=eq_G1,
        eq2=eq_G2,
        ts1=ts_G1,
        ts2=ts_G2,
        sol1=sol_G1,
        sol2=sol_G2,
        sol1_n = subset(sol_G1, sol_G1[,4]<=input$zd_1),
        sol1_d = subset(sol_G1, sol_G1[,4]>input$zd_1),
        sol2_n = subset(sol_G2, sol_G2[,4]<=input$zd_2),
        sol2_d = subset(sol_G2, sol_G2[,4]>input$zd_2),
        pdt=pdt,
        hue=hue
      )
    )
  })
  
  ###generate stoch sims
  simus_stoch_G1 <- reactive({
    par_G1 <- as.list(c(delta=input$delta_1,
                        phi=input$phi_1,
                        gamma=input$gamma_1,
                        alpha=input$alpha_1,
                        beta=input$beta_1,
                        u=input$u_1, 
                        v=input$v_1,
                        w=input$w_1, #this is now called k
                        psi=input$psi_1,
                        l=input$l_1,
                        omega=input$omega_1,
                        eta=input$eta_1,
                        ksi=input$ksi_1,
                        zcrit=input$zd_1 #also called zd
    )
    )
    
    eq_G1 <- equilibrium(par_G1)
    ###stoch sim
    n1 <- input$npts1
    splash_G1 <- c()
    pdt_stoch <- seq(0,input$tmax,by=20)
    
    par_ini_stoch_1=as.list(c(mu_x=10^input$mu_x1,
                              mu_y=input$mu_y1,
                              mu_z=input$mu_z1,
                              sigma_logx=input$var_x1,
                              sigma_logy=input$var_y1,
                              sigma_logz=input$var_z1))
    for(i in pdt_stoch){
      dtemp_G1 <- survsim(n1,par_G1,par_ini_stoch_1,tmax = i)
      #dtemp_G1 <- subset(dtemp_G1, dead==0)
      splash_G1 <- rbind(splash_G1,dtemp_G1)
      rm(dtemp_G1)
    }
    splash_dead_G1 <- subset(splash_G1, dead==1)
    splash_G1 <- subset(splash_G1, dead==0)
    ###setting up the colour object
    if(input$clrb){
      hue<-c("#2F4562","#C14C32")
    }else{
      hue <- c("#1D7874","#C39059")
      
    }
    
    obj_lst <- return(
      list(
        par_1=par_G1,
        eq1=eq_G1,
        stoch1 = splash_G1,
        stoch_dead1 = splash_dead_G1,
        pdt_stoch=pdt_stoch,
        hue=hue
      )
    )
    
  })
  simus_stoch_G2 <- reactive({
    par_G2 <- as.list(c(delta=input$delta_2,
                        phi=input$phi_2,
                        gamma=input$gamma_2,
                        alpha=input$alpha_2,
                        beta=input$beta_2,
                        u=input$u_2, 
                        v=input$v_2,
                        w=input$w_2, #this is now called k
                        psi=input$psi_2,
                        l=input$l_2,
                        omega=input$omega_2,
                        eta=input$eta_2,
                        ksi=input$ksi_2,
                        zcrit=input$zd_2 #also called zd
    )
    )
    
    eq_G2 <- equilibrium(par_G2)
    
    n2 <- input$npts2
    splash_G2 <- c()
    pdt_stoch <- seq(0,input$tmax,by=20)
    
    par_ini_stoch_2=as.list(c(mu_x=10^input$mu_x2,
                              mu_y=input$mu_y2,
                              mu_z=input$mu_z2,
                              sigma_logx=input$var_x2,
                              sigma_logy=input$var_y2,
                              sigma_logz=input$var_z2))
    
    for(i in pdt_stoch){
      dtemp_G2 <- survsim(n2,par_G2,par_ini_stoch_2,tmax = i)
      #dtemp_G2 <- subset(dtemp_G2, dead==0)
      splash_G2 <- rbind(splash_G2,dtemp_G2)
      rm(dtemp_G2)
    }
    splash_dead_G2 <- subset(splash_G2, dead==1)
    splash_G2 <- subset(splash_G2, dead==0)
    ###setting up the colour object
    if(input$clrb){
      hue<-c("#2F4562","#C14C32")
    }else{
      hue <- c("#1D7874","#C39059")#,"#F7441B"
    }
    
    obj_lst <- return(
      list(
        par_2=par_G2,
        eq2=eq_G2,
        stoch2 = splash_G2,
        stoch_dead2 = splash_dead_G2,
        pdt_stoch=pdt_stoch,
        hue=hue
      )
    )
  })
  simus_surv <- reactive({
    par_G1 <- as.list(c(delta=input$delta_1,
                        phi=input$phi_1,
                        gamma=input$gamma_1,
                        alpha=input$alpha_1,
                        beta=input$beta_1,
                        u=input$u_1, 
                        v=input$v_1,
                        w=input$w_1, #this is now called k
                        psi=input$psi_1,
                        l=input$l_1,
                        omega=input$omega_1,
                        eta=input$eta_1,
                        ksi=input$ksi_1,
                        zcrit=input$zd_1 #also called zd
    )
    )
    par_G2 <- as.list(c(delta=input$delta_2,
                        phi=input$phi_2,
                        gamma=input$gamma_2,
                        alpha=input$alpha_2,
                        beta=input$beta_2,
                        u=input$u_2, 
                        v=input$v_2,
                        w=input$w_2, #this is now called k
                        psi=input$psi_2,
                        l=input$l_2,
                        omega=input$omega_2,
                        eta=input$eta_2,
                        ksi=input$ksi_2,
                        zcrit=input$zd_2 #also called zd
    )
    )
    eq_G1 <- equilibrium(par_G1)
    eq_G2 <- equilibrium(par_G2)
    
    n <- input$nb_surv
    splash_G1 <- c()
    splash_G2 <- c()
    #pdt_stoch <- seq(0,input$tmax,by=20)
    
    par_ini_surv1=as.list(c(mu_x=10^input$mu_x1,
                            mu_y=input$mu_y1,
                            mu_z=input$mu_z1,
                            sigma_logx=input$var_x1,
                            sigma_logy=input$var_y1,
                            sigma_logz=input$var_z1))
    
    par_ini_surv2=as.list(c(mu_x=10^input$mu_x2,
                            mu_y=input$mu_y2,
                            mu_z=input$mu_z2,
                            sigma_logx=input$var_x2,
                            sigma_logy=input$var_y2,
                            sigma_logz=input$var_z2))
    
    surv_G1 <- survsim(n,par_G1,par_ini_surv1,tmax = input$tmax)
    surv_G2 <- survsim(n,par_G2,par_ini_surv2,tmax = input$tmax)
    
    surv_dead_G1 <- subset(surv_G1, dead==1)
    surv_dead_G2 <- subset(surv_G2, dead==1)
    ###setting up the colour object
    if(input$clrb){
      hue<-c("#2F4562","#C14C32")
    }else{
      hue <- c("#1D7874","#C39059")#,"#F7441B"
    }
    
    obj_lst <- return(
      list(
        par_1=par_G2,
        par_2=par_G2,
        eq1=eq_G1,
        eq2=eq_G2,
        surv1 = surv_G1,
        surv2 = surv_G2,
        surv_dead1 = surv_dead_G1,
        surv_dead2 = surv_dead_G2,
        hue=hue
      )
    )
  })
  ###generate table of difference between par sets
  output$var_pars_panel1 <- renderUI({
    validate(
      need(input$delta_1 != "", "Loading...") # display custom message in need
    )
    par_1 <- as.list(c(delta=input$delta_1,
                       phi=input$phi_1,
                       gamma=input$gamma_1,
                       alpha=input$alpha_1,
                       beta=input$beta_1,
                       u=input$u_1, 
                       v=input$v_1,
                       w=input$w_1, #this is now called k
                       psi=input$psi_1,
                       l=input$l_1,
                       omega=input$omega_1,
                       eta=input$eta_1,
                       ksi=input$ksi_1,
                       zcrit=input$zd_1
    )
    )
    par_2 <- as.list(c(delta=input$delta_2,
                       phi=input$phi_2,
                       gamma=input$gamma_2,
                       alpha=input$alpha_2,
                       beta=input$beta_2,
                       u=input$u_2, 
                       v=input$v_2,
                       w=input$w_2, #this is now called k
                       psi=input$psi_2,
                       l=input$l_2,
                       omega=input$omega_2,
                       eta=input$eta_2,
                       ksi=input$ksi_2,
                       zcrit=input$zd_2
    )
    )
    
    Variation <- unlist(par_1)-unlist(par_2)
    delta_par <- as.data.frame(Variation)
    row.names(delta_par) <- gsub("ksi","xi",rownames(delta_par))
    row.names(delta_par) <- paste("\\",row.names(delta_par),sep="")
    row.names(delta_par) <- gsub("\\\\(.)$","\\1",rownames(delta_par))
    row.names(delta_par) <- gsub("\\\\(zcrit)$","zd",rownames(delta_par))
    
    ##last object returned
    
    delta_par <- subset(delta_par, Variation!=0)
    colnames(delta_par) <- "\\Delta (G.1 - G.2)"
    
    LaTeXtab <- print(xtable(delta_par, align=rep("c", ncol(delta_par)+1)), 
                      floating=FALSE, tabular.environment="array", comment=FALSE, 
                      print.results=FALSE, 
                      sanitize.rownames.function = function(x) x,
                      sanitize.colnames.function = function(x) x)
    tagList(
      withMathJax(),
      HTML(paste0("$$", LaTeXtab, "$$"))
    )
    
  })
  output$var_pars_panel2 <- renderUI({
    validate(
      need(input$gamma_G1_chckbox != "", "Please go to the settings tab, parameters need to be initialized to compute the summary data.") # display custom message in need
    )
    
    validate(
      need(input$delta_1 != "", "Loading...") # display custom message in need
    )
    par_1 <- as.list(c(delta=input$delta_1,
                       phi=input$phi_1,
                       gamma=input$gamma_1,
                       alpha=input$alpha_1,
                       beta=input$beta_1,
                       u=input$u_1, 
                       v=input$v_1,
                       w=input$w_1, #this is now called k
                       psi=input$psi_1,
                       l=input$l_1,
                       omega=input$omega_1,
                       eta=input$eta_1,
                       ksi=input$ksi_1,
                       zcrit=input$zd_1
    )
    )
    par_2 <- as.list(c(delta=input$delta_2,
                       phi=input$phi_2,
                       gamma=input$gamma_2,
                       alpha=input$alpha_2,
                       beta=input$beta_2,
                       u=input$u_2, 
                       v=input$v_2,
                       w=input$w_2, #this is now called k
                       psi=input$psi_2,
                       l=input$l_2,
                       omega=input$omega_2,
                       eta=input$eta_2,
                       ksi=input$ksi_2,
                       zcrit=input$zd_2
    )
    )
    
    Variation <- unlist(par_1)-unlist(par_2)
    delta_par <- as.data.frame(Variation)
    row.names(delta_par) <- gsub("ksi","xi",rownames(delta_par))
    row.names(delta_par) <- paste("\\",row.names(delta_par),sep="")
    row.names(delta_par) <- gsub("\\\\(.)$","\\1",rownames(delta_par))
    row.names(delta_par) <- gsub("\\\\(zcrit)$","zd",rownames(delta_par))
    ##last object returned
    names_par <-c("Pathogen\\, destruction\\, rate\\,",
                  "Defense\\, decay\\ rate\\,",
                  "Constitutive\\, immunity\\,",
                  "Immune\\, activation\\,",
                  "Immune\\, saturation\\,",
                  "Activation\\, shape\\, parameter\\,",
                  "Saturation\\, shape\\, parameter\\,",
                  "Down-regulation\\, shape\\, parameter\\,",
                  "Strength\\, of\\, damage\\, effect\\,",
                  "Damage\\, effect\\, shape\\, parameter\\,",
                  "Damage\\, production\\, rate\\, from\\, pathogens\\,",
                  "Damage\\, production\\, rate\\, from\\, defense\\,",
                  "Damage\\, repair\\, rate\\,",
                  "Damage\\, threshold\\,"
    )
    
    for(i in 1:length(row.names(delta_par))){
      row.names(delta_par)[i] <- paste0(names_par[i]," (",row.names(delta_par)[i],")")
    }
    
    delta_par <- subset(delta_par, Variation!=0)
    colnames(delta_par) <- "\\Delta (G.1 - G.2)"
    
    LaTeXtab <- print(xtable(delta_par, align=rep("c", ncol(delta_par)+1)), 
                      floating=FALSE, tabular.environment="array", comment=FALSE, 
                      print.results=FALSE, 
                      sanitize.rownames.function = function(x) x,
                      sanitize.colnames.function = function(x) x)
    tagList(
      withMathJax(),
      HTML(paste0("$$", LaTeXtab, "$$"))
    )
    
  })#
  output$var_pars_panel3 <- renderUI({
    validate(
      need(input$gamma_G1_chckbox != "", "Please go to the settings tab, parameters need to be initialized to compute the summary data.") # display custom message in need
    )
    
    validate(
      need(input$delta_1 != "", "Loading...") # display custom message in need
    )
    par_1 <- as.list(c(delta=input$delta_1,
                       phi=input$phi_1,
                       gamma=input$gamma_1,
                       alpha=input$alpha_1,
                       beta=input$beta_1,
                       u=input$u_1, 
                       v=input$v_1,
                       w=input$w_1, #this is now called k
                       psi=input$psi_1,
                       l=input$l_1,
                       omega=input$omega_1,
                       eta=input$eta_1,
                       ksi=input$ksi_1,
                       zcrit=input$zd_1
    )
    )
    par_2 <- as.list(c(delta=input$delta_2,
                       phi=input$phi_2,
                       gamma=input$gamma_2,
                       alpha=input$alpha_2,
                       beta=input$beta_2,
                       u=input$u_2, 
                       v=input$v_2,
                       w=input$w_2, #this is now called k
                       psi=input$psi_2,
                       l=input$l_2,
                       omega=input$omega_2,
                       eta=input$eta_2,
                       ksi=input$ksi_2,
                       zcrit=input$zd_2
    )
    )
    
    Variation <- unlist(par_1)-unlist(par_2)
    delta_par <- as.data.frame(Variation)
    row.names(delta_par) <- gsub("ksi","xi",rownames(delta_par))
    row.names(delta_par) <- paste("\\",row.names(delta_par),sep="")
    row.names(delta_par) <- gsub("\\\\(.)$","\\1",rownames(delta_par))
    row.names(delta_par) <- gsub("\\\\(zcrit)$","zd",rownames(delta_par))
    ##last object returned
    names_par <-c("Pathogen\\, destruction\\, rate\\,",
                  "Defense\\, decay\\ rate\\,",
                  "Constitutive\\, immunity\\,",
                  "Immune\\, activation\\,",
                  "Immune\\, saturation\\,",
                  "Activation\\, shape\\, parameter\\,",
                  "Saturation\\, shape\\, parameter\\,",
                  "Down-regulation\\, shape\\, parameter\\,",
                  "Strength\\, of\\, damage\\, effect\\,",
                  "Damage\\, effect\\, shape\\, parameter\\,",
                  "Damage\\, production\\, rate\\, from\\, pathogens\\,",
                  "Damage\\, production\\, rate\\, from\\, defense\\,",
                  "Damage\\, repair\\, rate\\,",
                  "Damage\\, threshold\\,"
    )
    
    for(i in 1:length(row.names(delta_par))){
      row.names(delta_par)[i] <- paste0(names_par[i]," (",row.names(delta_par)[i],")")
    }
    
    delta_par <- subset(delta_par, Variation!=0)
    colnames(delta_par) <- "\\Delta (G.1 - G.2)"
    
    LaTeXtab <- print(xtable(delta_par, align=rep("c", ncol(delta_par)+1)), 
                      floating=FALSE, tabular.environment="array", comment=FALSE, 
                      print.results=FALSE, 
                      sanitize.rownames.function = function(x) x,
                      sanitize.colnames.function = function(x) x)
    tagList(
      withMathJax(),
      HTML(paste0("$$", LaTeXtab, "$$"))
    )
    
  })
  
  ###make plots
  observe({
    output$dyn_plotly_x <- renderPlotly({
      validate(
        need(input$gamma_G1_chckbox != "", "Please go to the settings tab first in order for the graphs to compute.") # display custom message in need
      )
      
      validate(
        need(input$x_ini1 != "", "Loading...") # display custom message in need
      )
      ##load sim data
      dyn_S <- simus()
      sol1 <- dyn_S$sol1
      sol1_n <- dyn_S$sol1_n
      sol1_d <- dyn_S$sol1_d
      
      sol2 <- dyn_S$sol2
      sol2_n <- dyn_S$sol2_n
      sol2_d <- dyn_S$sol2_d
      
      ts1 <- dyn_S$ts1
      ts2 <- dyn_S$ts2
      
      hue <- dyn_S$hue
      ##plot
      f_X <- ggplot() +
        geom_line(sol1_n, mapping=aes(x=t,y=x, group=1,
                                      text = sprintf("Genotype 1<br>Time: %s<br>P. load (<i>x</i>): %s",
                                                     t, scales::scientific(x,digits = 3,scale=1,trim=TRUE)
                                      )
        ),size=1,colour=hue[1]) +
        geom_line(sol2_n, mapping=aes(x=t,y=x, group = 2,
                                      text = sprintf("Genotype 2<br>Time: %s<br>P. load (<i>x</i>): %s",
                                                     t, scales::scientific(x,digits = 3,scale=1,trim=TRUE)
                                      )
        ),size=1,colour=hue[2]) +
        labs(x="Time (<i>t</i>)",y="Pathogen load (<i>x</i>)") +
        #coord_cartesian(xlim = range(sol1[,1]),expand = FALSE) +
        theme(axis.title = element_text(size=14), ##changes axis names sizes
              axis.text = element_text(size=12) ## changes axis marks sizes
        )
      if(nrow(sol1_d)!=0){
        f_X <- f_X + geom_line(sol1_d, mapping=aes(t,x, group = 3,
                                                   text = sprintf("G.1 after death<br>Time: %s<br>P. load (<i>x</i>): %s",
                                                                  t, scales::scientific(x,digits = 3,scale=1,trim=TRUE)
                                                   )
        ),linetype="dashed",colour=hue[1]) +
          geom_point(mapping=aes(ts1[4],ts1[1],#group =3,
                                 text = sprintf("Death of G.1<br>Time to death: %s<br>PLUD: %s",
                                                round(ts1[4],digits = 2), scales::scientific(ts1[1],digits = 3,scale=1,trim=TRUE)
                                 )
          ),colour=hue[1])
      }
      if(nrow(sol2_d)!=0){
        f_X <- f_X + geom_line(sol2_d, mapping=aes(t,x, group = 3,
                                                   text = sprintf("G.2 after death<br>Time: %s<br>P. load (<i>x</i>): %s",
                                                                  t, scales::scientific(x,digits = 3,scale=1,trim=TRUE)
                                                   )
        ),linetype="dashed",colour=hue[2]) +
          geom_point(mapping=aes(ts2[4],ts2[1],                   # ,group =3,
                                 text = sprintf("Death of G.2<br>Time to death: %s<br>PLUD: %s",
                                                round(ts2[4],digit=2), scales::scientific(ts2[1],digits = 3,scale=1,trim=TRUE)
                                 )
          ),colour=hue[2])
      }
      if(input$scale_dynx == "log10"){
        f_X <- f_X + scale_x_log10()
      }
      if(input$scale_dyny == "log2"){
        f_X <- f_X + scale_y_continuous(trans = 'log2', breaks = c(1e-6,.1,1))
      }
      if(input$scale_dyny == "log10"){
        f_X <- f_X + scale_y_log10()
      }  
      
      f_X <- ggplotly(f_X, tooltip="text") %>% config(displayModeBar = FALSE) #%>% layout(hovermode='x')#%>% style(hovertext=list(digits=1))
      f_X
    })
    
    output$dyn_plotly_y <- renderPlotly({
      validate(
        need(input$gamma_G1_chckbox != "", "") # display custom message in need
      )
      
      validate(
        need(input$x_ini1 != "", "") # display custom message in need
      )
      ##load sim data
      dyn_S <- simus()
      sol1 <- dyn_S$sol1
      sol1_n <- dyn_S$sol1_n
      sol1_d <- dyn_S$sol1_d
      
      sol2 <- dyn_S$sol2
      sol2_n <- dyn_S$sol2_n
      sol2_d <- dyn_S$sol2_d
      
      hue <- dyn_S$hue
      ##plot
      f_Y <- ggplot() +
        geom_line(sol1_n, mapping=aes(x=t,y=y, group=1,
                                      text = sprintf("Genotype 1<br>Time: %s<br>Defense (<i>y</i>): %s",
                                                     t, scales::scientific(y,digits = 3,scale=1,trim=TRUE)
                                      )
        ),size=1,colour=hue[1]) +
        geom_line(sol2_n, mapping=aes(t,y, group = 2,
                                      text = sprintf("Genotype 2<br>Time: %s<br>Defense (<i>y</i>): %s",
                                                     t, scales::scientific(y,digits = 3,scale=1,trim=TRUE)
                                      )
        ),size=1,colour=hue[2]) +
        labs(x="Time (<i>t</i>)",y="Host defense (<i>y</i>)") +
        theme(axis.title = element_text(size=14), ##changes axis names sizes
              axis.text = element_text(size=12) ## changes axis marks sizes
        )
      if(nrow(sol1_d)!=0){
        f_Y <- f_Y + geom_line(sol1_d, mapping=aes(t,y, group = 3,
                                                   text = sprintf("G.1 after death<br>Time: %s<br>Defense (<i>y</i>): %s",
                                                                  t, scales::scientific(y,digits = 3,scale=1,trim=TRUE)
                                                   )
        ),linetype="dashed",colour=hue[1]) 
      }
      if(nrow(sol2_d)!=0){
        f_Y <- f_Y + geom_line(sol2_d, mapping=aes(t,y, group = 3,
                                                   text = sprintf("G.2 after death<br>Time: %s<br>Defense (<i>y</i>): %s",
                                                                  t, scales::scientific(y,digits = 3,scale=1,trim=TRUE)
                                                   )
        ),linetype="dashed",colour=hue[2])
      }
      if(input$scale_dynx == "log10"){
        f_Y <- f_Y + scale_x_log10()
      }
      if(input$scale_dyny == "log2"){
        f_Y <- f_Y + scale_y_continuous(trans = 'log2', breaks = c(1e-6,.1,1))
      }
      if(input$scale_dyny == "log10"){
        f_Y <- f_Y + scale_y_log10()
      }  
      f_Y <- ggplotly(f_Y, tooltip="text") %>% config(displayModeBar = FALSE)
      f_Y
      
    })
    
    output$dyn_plotly_z <- renderPlotly({
      validate(
        need(input$gamma_G1_chckbox != "", "") # display custom message in need
      )
      
      validate(
        need(input$x_ini1 != "", "") # display custom message in need
      )
      ##load sim data
      dyn_S <- simus()
      sol1 <- dyn_S$sol1
      sol1_n <- dyn_S$sol1_n
      sol1_d <- dyn_S$sol1_d
      
      sol2 <- dyn_S$sol2
      sol2_n <- dyn_S$sol2_n
      sol2_d <- dyn_S$sol2_d
      
      hue <- dyn_S$hue
      
      ##plot
      f_Z <- ggplot() +
        geom_line(sol1_n, mapping=aes(x=t,y=z, group=1,
                                      text = sprintf("Genotype 1<br>Time: %s<br>Damage (<i>z</i>): %s",
                                                     t, scales::scientific(z,digits = 3,scale=1,trim=TRUE)
                                      )
        ),size=1,colour=hue[1]) +
        geom_line(sol2_n, mapping=aes(t,z, group = 2,
                                      text = sprintf("Genotype 2<br>Time: %s<br>Damage (<i>z</i>): %s",
                                                     t, scales::scientific(z,digits = 3,scale=1,trim=TRUE)
                                      )
        ),size=1,colour=hue[2]) +
        labs(x="Time (<i>t</i>)",y="Host damage (<i>z</i>)") +
        theme(axis.title = element_text(size=14), ##changes axis names sizes
              axis.text = element_text(size=12) ## changes axis marks sizes
        )
      if(nrow(sol1_d)!=0){
        f_Z <- f_Z + geom_line(sol1_d, mapping=aes(t,z, group = 3,
                                                   text = sprintf("G.1 after death<br>Time: %s<br>Damage (<i>z</i>): %s",
                                                                  t, scales::scientific(z,digits = 3,scale=1,trim=TRUE)
                                                   )
        ),linetype="dashed",colour=hue[1])
      }
      if(nrow(sol2_d)!=0){
        f_Z <- f_Z + geom_line(sol2_d, mapping=aes(t,z, group = 3,
                                                   text = sprintf("G.2 after death<br>Time: %s<br>Damage (<i>z</i>): %s",
                                                                  t, scales::scientific(z,digits = 3,scale=1,trim=TRUE)
                                                   )
        ),linetype="dashed",colour=hue[2])
      }
      if(input$zd_1 == input$zd_2){
        f_Z <- f_Z + geom_hline(aes(yintercept = input$zd_1,text=paste("<i>zcrit</i>")),linetype="dotted",size=.3)
      }else{
        f_Z <- f_Z + geom_hline(aes(yintercept = input$zd_1,text=paste("<i>zcrit (G.1)</i>")),linetype="dotted",size=.3) +
          geom_hline(aes(yintercept = input$zd_2,text=paste("<i>zcrit (G.2)</i>")),linetype="dotted",size=.3)
      }
      if(input$scale_dynx == "log10"){
        f_Z <- f_Z + scale_x_log10()
      }
      if(input$scale_dyny == "log2"){
        f_Z <- f_Z + scale_y_continuous(trans = 'log2', breaks = c(1e-6,.1,1))
      }
      if(input$scale_dyny == "log10"){
        f_Z <- f_Z + scale_y_log10()
      }  
      f_Z <- ggplotly(f_Z, tooltip="text") %>% config(displayModeBar = FALSE)
      f_Z
    })
  })
  output$phase_plot_3D <- renderPlotly({
    validate(
      need(input$gamma_G1_chckbox != "", "Please go to the settings tab first in order for the graphs to compute.") # display custom message in need
    )
    
    validate(
      need(input$delta_1 != "", "Loading...") # display custom message in need
    )
    
    dyn_S <- simus()
    sol1 <- dyn_S$sol1
    sol1_n <- dyn_S$sol1_n
    sol1_d <- dyn_S$sol1_d
    eq1 <- dyn_S$eq1
    
    sol2 <- dyn_S$sol2
    sol2_n <- dyn_S$sol2_n
    sol2_d <- dyn_S$sol2_d
    eq2 <- dyn_S$eq2
    
    hue <- dyn_S$hue
    
    fig3D <- plot_ly()
    
    fig3D <- fig3D %>% add_trace(data=sol1_n,
                                 x=sol1_n[,2],
                                 y=sol1_n[,3],
                                 z=sol1_n[,4],
                                 name="Traj. G1",
                                 type="scatter3d",
                                 mode="lines",
                                 line=list(width=6,color=hue[1]) #default size appears to be 2
    )
    
    fig3D <- fig3D %>% add_trace(data=sol1_n,
                                 x=sol1_d[,2],
                                 y=sol1_d[,3],
                                 z=sol1_d[,4],
                                 name="Traj. G1 after death",
                                 type="scatter3d",
                                 mode="lines",
                                 line=list(width=3,color=hue[1],dash="dash") #default size appears to be 2
    )
    
    
    fig3D <- fig3D %>% add_trace(data=sol2_n,
                                 x=sol2_n[,2],
                                 y=sol2_n[,3],
                                 z=sol2_n[,4],
                                 name="Traj. G2",
                                 mode="lines",
                                 line=list(width=6, color=hue[2])
    )
    
    fig3D <- fig3D %>% add_trace(data=sol2_d,
                                 x=sol2_d[,2],
                                 y=sol2_d[,3],
                                 z=sol2_d[,4],
                                 name="Traj. G2 after death",
                                 type="scatter3d",
                                 mode="lines",
                                 line=list(width=3,color=hue[2],dash="dash") #default size appears to be 2
    )
    ## adding ed values for genotype 1
    fig3D <- fig3D %>% add_trace(data=eq1[-1,],
                                 x=~x,
                                 y=~y,
                                 z=~z,
                                 name="equilibria G1",
                                 mode="markers",
                                 marker = list(size=2,color = hue[1]))
    
    # adding ed values for genotype 2
    fig3D <- fig3D %>% add_trace(data=eq2[-1,],
                                 x=~x,
                                 y=~y,
                                 z=~z,
                                 name="equilibria G2",
                                 mode="markers",
                                 marker = list(size=2,color = hue[2]))
    
    fig3D <- fig3D %>% layout(paper_bgcolor='#222222',
                              scene=list(
                                xaxis = list(title="Pathogen load (x)",
                                             color="#fafafa",
                                             showline=FALSE,
                                             showticklabels=FALSE,
                                             showgrid=TRUE,
                                             autorange=TRUE),
                                yaxis = list(title="Host defense (y)",
                                             color="#fafafa",
                                             showline=FALSE,
                                             showticklabels=FALSE,
                                             showgrid=TRUE,
                                             autorange=TRUE),
                                zaxis = list(title="Host damage (z)",
                                             color="#fafafa",
                                             showline=FALSE,
                                             showticklabels=FALSE,
                                             showgrid=TRUE,
                                             autorange=TRUE)
                              ))
  })
  output$dyn_plot_3D <- renderPlotly({
    dyn_S <- simus()
    sol1 <- dyn_S$sol1
    sol1_n <- dyn_S$sol1_n
    sol1_d <- dyn_S$sol1_d
    eq1 <- dyn_S$eq1
    
    sol2 <- dyn_S$sol2
    sol2_n <- dyn_S$sol2_n
    sol2_d <- dyn_S$sol2_d
    eq2 <- dyn_S$eq2
    
    fig3D <- plot_ly(data=sol1,
                     x=sol1[,2],
                     y=sol1[,3],
                     z=sol1[,4],
                     name="Traj. G1",
                     type="scatter3d",
                     mode="lines"
    )
    
    fig3D <- fig3D %>% add_trace(data=sol2,
                                 x=sol2[,2],
                                 y=sol2[,3],
                                 z=sol2[,4],
                                 name="Traj. G2",
                                 mode="lines"
    )
    ## adding ed values for genotype 1
    fig3D <- fig3D %>% add_trace(data=eq1[-1,],
                                 x=~x,
                                 y=~y,
                                 z=~z,
                                 name="equilibria G1",
                                 mode="markers",
                                 marker = list(size=2,color = "steelblue"))
    
    # adding ed values for genotype 2
    fig3D <- fig3D %>% add_trace(data=eq2[-1,],
                                 x=~x,
                                 y=~y,
                                 z=~z,
                                 name="equilibria G2",
                                 mode="markers",
                                 marker = list(size=2,color = "orange"))
    
    fig3D <- fig3D %>% layout(paper_bgcolor='#222222',
                              scene=list(
                                xaxis = list(title="Pathogen load (x)",
                                             color="#fafafa",
                                             showline=FALSE,
                                             showticklabels=FALSE,
                                             showgrid=TRUE,
                                             autorange=TRUE),
                                yaxis = list(title="Level of defense (y)",
                                             color="#fafafa",
                                             showline=FALSE,
                                             showticklabels=FALSE,
                                             showgrid=TRUE,
                                             autorange=TRUE),
                                zaxis = list(title="Amount of damage (z)",
                                             color="#fafafa",
                                             showline=FALSE,
                                             showticklabels=FALSE,
                                             showgrid=TRUE,
                                             autorange=TRUE)
                              ))
  })
  output$stoch_plotly_G1 <- renderPlotly({
    validate(
      need(input$gamma_G1_chckbox != "", "Please go to the settings tab first in order for the graphs to compute.") # display custom message in need
    )
    
    validate(
      need(input$mu_x1 != "", "Loading...") # display custom message in need
    )
    ##load sim data
    stc_S <- simus_stoch_G1()
    stoch1 <- stc_S$stoch1
    #stoch2 <- stc_S$stoch2
    
    hue <- stc_S$hue
    ##plot
    f_sX <- ggplot() +
      geom_point(stoch1, mapping=aes(x=time,y=x,group=1,
                                     text = sprintf("Genotype 1<br>Time: %s<br>P. load (<i>x</i>): %s",
                                                    time, scales::scientific(x,digits = 3,scale=1,trim=TRUE)
                                     )
      ),size=1,colour=hue[1]) +
      labs(x="Time (<i>t</i>)",y="Pathogen load (<i>x</i>)") +
      #coord_cartesian(xlim = range(sol1[,1]),expand = FALSE) +
      theme(axis.title = element_text(size=14), ##changes axis names sizes
            axis.text = element_text(size=12) ## changes axis marks sizes
      )
    if(input$scale_stoch1x == "log10"){
      f_sX <- f_sX + scale_x_log10()
    }
    if(input$scale_stoch1y == "log2"){
      f_sX <- f_sX + scale_y_continuous(trans = 'log2', breaks = c(1e-6,.1,1))
    }
    if(input$scale_stoch1y == "log10"){
      f_sX <- f_sX + scale_y_log10()
    }  
    f_sX <- ggplotly(f_sX, tooltip="text") %>% config(displayModeBar = FALSE) #%>% layout(hovermode='x')#%>% style(hovertext=list(digits=1))
    f_sX
  })
  
  output$stoch_plotly_G2 <- renderPlotly({
    validate(
      need(input$gamma_G2_chckbox != "", "") # display custom message in need
    )
    
    validate(
      need(input$mu_x2 != "", "Loading...") # display custom message in need
    )
    ##load sim data
    stc_S <- simus_stoch_G2()
    #stoch1 <- stc_S$stoch1
    stoch2 <- stc_S$stoch2
    
    hue <- stc_S$hue
    ##plot
    f_sX <- ggplot() +
      geom_point(stoch2, mapping=aes(x=time,y=x,group=1,
                                     text = sprintf("Genotype 2<br>Time: %s<br>P. load (<i>x</i>): %s",
                                                    time, scales::scientific(x,digits = 3,scale=1,trim=TRUE)
                                     )
      ),size=1,colour=hue[2]) +
      labs(x="Time (<i>t</i>)",y="Pathogen load (<i>x</i>)") +
      #coord_cartesian(xlim = range(sol1[,1]),expand = FALSE) +
      theme(axis.title = element_text(size=14), ##changes axis names sizes
            axis.text = element_text(size=12) ## changes axis marks sizes
      )
    if(input$scale_stoch2x == "log10"){
      f_sX <- f_sX + scale_x_log10()
    }
    if(input$scale_stoch2y == "log2"){
      f_sX <- f_sX + scale_y_continuous(trans = 'log2', breaks = c(1e-6,.1,1))
    }
    if(input$scale_stoch2y == "log10"){
      f_sX <- f_sX + scale_y_log10()
    }  
    f_sX <- ggplotly(f_sX, tooltip="text") %>% config(displayModeBar = FALSE) #%>% layout(hovermode='x')#%>% style(hovertext=list(digits=1))
    f_sX
  })
  
  # output$hist_dead <- renderPlotly({
  #   validate(need(input$gamma_G2_chckbox != "", "") # display custom message in need)
  #            
  #            validate(need(input$mu_x2 != "", "Loading...") # display custom message in need)
  #                     ##load sim data
  #                     stc_S1 <- simus_stoch_G1()
  #                     stc_S2 <- simus_stoch_G2()
  #                     #stoch1 <- stc_S$stoch1
  #                     stoch_dead1 <- stc_S1$stoch_dead1
  #                     stoch_dead2 <- stc_S2$stoch_dead2
  #                     
  #                     hue <- stc_S1$hue
  #                     
  #                     H <- ggplot() +
  #                       geom_histogram(
  #                         stoch_dead1,
  #                         mapping = aes(x = time, y = ..count.., label = ..count..),
  #                         fill = hue[1],
  #                         binwidth = 20
  #                       ) +
  #                       geom_histogram(
  #                         stoch_dead2,
  #                         mapping = aes(
  #                           x = time,
  #                           y = -..count..,
  #                           label = ..count..
  #                         ),
  #                         fill = hue[2],
  #                         binwidth = 20
  #                       ) +
  #                       scale_y_continuous(breaks = seq(-300, 300, by = 25),
  #                                          labels = abs(seq(-300, 300, by = 25))) +
  #                       labs(x = "Time (<i>t</i>)", y = "Dead hosts") +
  #                       #coord_cartesian(xlim = range(sol1[,1]),expand = FALSE) +
  #                       theme(axis.title = element_text(size = 14),
  #                             ##changes axis names sizes
  #                             axis.text = element_text(size = 12) ## changes axis marks sizes)
  #                             H <-
  #                               ggplotly(H, tooltip = 'label') %>% config(displayModeBar = FALSE) #%>% layout(hovermode='x')#%>% style(hovertext=list(digits=1))
  #                             H
  # })
  
  output$surv_dead <- renderPlotly({
    validate(
      need(input$gamma_G2_chckbox != "", "") # display custom message in need
    )
    
    validate(
      need(input$mu_x2 != "", "Loading...") # display custom message in need
    )
    ##load sim data
    surv_d <- simus_surv()
    #stoch1 <- stc_S$stoch1
    surv_dead1 <- surv_d$surv_dead1
    surv_dead2 <- surv_d$surv_dead2
    
    hue <- surv_d$hue
    
    ##reorder
    
    S <- ggplot()
    
    if(nrow(surv_dead1)>0){
      surv_dead1 <- surv_dead1[order(surv_dead1$time),]
      surv_dead1$surv_density <- NA
      for(i in 1:nrow(surv_dead1)){
        surv_dead1$surv_density[i] <- 1 - (i/input$nb_surv)
      }
      surv_dead1 <- rbind(c(rep(0,8),1),surv_dead1)
      S <- S + geom_step(surv_dead1, mapping = aes(x=time,y=surv_density,group=1,
                                                   text = sprintf("Genotype 1<br>Time: %s<br>Proportion: %s",
                                                                  round(time,3), round(surv_density,3)
                                                   )
      ),colour=hue[1])
      
    }else{
      S <- S + geom_hline(mapping = aes(yintercept=1,text="No death for G.1"),colour=hue[1])
    }
    
    if(nrow(surv_dead2)>0){
      surv_dead2 <- surv_dead2[order(surv_dead2$time),]
      surv_dead2$surv_density <- NA
      for(i in 1:nrow(surv_dead2)){
        surv_dead2$surv_density[i] <- 1 - (i/input$nb_surv) 
      }
      surv_dead2 <- rbind(c(rep(0,8),1),surv_dead2)
      S <- S + geom_step(surv_dead2, mapping = aes(x=time,y=surv_density,group=2,
                                                   text = sprintf("Genotype 2<br>Time: %s<br>Proportion: %s",
                                                                  round(time,3), round(surv_density,3)
                                                   )
      ),colour=hue[2])
    }else{
      S <- S + geom_hline(mapping = aes(yintercept=1,text="No death for G.2"),colour=hue[2])
    }
    S <- S +labs(x="Time (<i>t</i>)",y="Proportion surviving") +
      #coord_cartesian(xlim = range(sol1[,1]),expand = FALSE) +
      theme(axis.title = element_text(size=14), ##changes axis names sizes
            axis.text = element_text(size=12) ## changes axis marks sizes
      )
    # if(input$scale_survx == "log10"){
    #   S <- S + scale_x_log10()
    # }
    # if(input$scale_survy == "log10"){
    #   S <- S + scale_y_log10()
    # }
    S <- ggplotly(S, tooltip='text') %>% config(displayModeBar = FALSE) #%>% layout(hovermode='x')#%>% style(hovertext=list(digits=1))
    S
    S
  })
  ###text outputs
  output$splash_dyn <- renderUI(HTML("<ul>
                                 <li>A dashed line signifies that the host is dead.</li>
                                 <li>When this happens, a point is plotted to represent the Pathogen Load Upon Death (PLUD).</li>
                                 </ul>
                                 "
  ))
  output$splash_stoch <- renderUI(HTML("<ul>
                                 <li>Here we sample the inoculum size around the mean to introduce tiny variations in the initial state of the system.</li>
                                 <li>This is used to mimic data we would get experimentally.</li>
                                 <li>Individuals that die are not represented on the scatter plots.</li>
                                 </ul>

                                 "
  ))
}

# Run the application 
shinyApp(ui = ui, server = server)
