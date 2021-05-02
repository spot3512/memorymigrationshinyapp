library(shiny)
ui <- fluidPage(
  
  h1("Memory Migration model"),
  fluidRow(
    column(3,
           h3("Migratory Population"),
           radioButtons(inputId = "world",
                        label = "Initial distribution of population in year 0", 
                        choices = c("Optimal" = "world_optimal", "Gaussian" = "world_gaussian", "Sinusoidal" = "world_sinusoidal"), inline = TRUE),
           numericInput(inputId = "years", label = "Duration of simulation:", value = 5, step = 1),
           numericInput(inputId = "threshold", label = "Threshold of Similarity", value = 0.999, min = 0, max = 1, step = 0.01),
           numericInput(inputId = "alphabeta",
                        label = "Total taxis",
                        value = 100, min = 0, step = 0.01),
           sliderInput(inputId = "kappa",
                       label = "Proportion resource following vs. memory",
                       value = .50, min = 0, max = 1),
           sliderInput(inputId = "gamma",
                       label = "Proportion reference memory",
                       value = 1, min = 0, max = 1),
           numericInput(inputId = "epsilon",
                        label = "Diffusion Parameter",
                        value = 1, min = 0, step = 0.01)
    ),
    column(9,
           h3("Migratory Population"),
           plotOutput("Image", height = "800px"),
           tableOutput("Indices")),
    column(3,
           h3("Resource"),
           sliderInput(inputId = "x.sd",
                       label = "Resource Space Distribution",
                       value = 3, min = 0, max = 15),
           sliderInput(inputId = "t.sd",
                       label = "Resource Time Distribution",
                       value = 3, min = 0, max = 15),
           numericInput(inputId = "mu.x0", label = "Initial Resource Position in Space", value = 80, step = 1),
           numericInput(inputId = "mu.t0", label = "Initial Resource Position in Time", value = 25, step = 1),
           numericInput(inputId = "beta.x", label = "Resource Change in Space", value = 0, step = 1),
           numericInput(inputId = "beta.t", label = "Resource Change in Time ", value = 0, step = 1),
           radioButtons(inputId = "resource",
                        label = "Type of resource", 
                        choices = c("Island" = "resources_island", "Drifting" = "resources_drifting"), inline = TRUE),
           
           actionButton(inputId = "run", 
                        label = "Run Model")),
    column(9,
           h3("Resource"),
           plotOutput("Resourceimage", height = "800px"))
    
  ))


server <- function(input, output) {
  pcks <- c("shiny","sf","ggplot2","magrittr","plyr", "gplots", "memorymigration")
  lapply(pcks, require, character = TRUE)
  
  simulation <- eventReactive(input$run, {
    if(input$world == "world_optimal"){
      world <- getOptimalPop(tau=100, X.min = 0, X.max = 100, dx=1, 
                             x.peak=as.numeric(input$mu.x0), t.peak=as.numeric(input$mu.t0), 
                             x.sd=as.numeric(input$x.sd), t.sd=as.numeric(input$t.sd))
    }
    else{
      data(world)
      world <- get(input$world)
    } 
    
    if(input$resource == "resources_island"){
      par0 <- getCCpars(mu_x0 = as.numeric(input$mu.x0), 
                        mu_t0 = as.numeric(input$mu.t0),
                        beta_x = as.numeric(input$beta.x),
                        beta_t = as.numeric(input$beta.t),
                        n.years = as.numeric(input$years),
                        sigma_x = as.numeric(input$x.sd),
                        sigma_t = as.numeric(input$t.sd))
      
      Resource.CC <- aaply(par0, 1, function(p) getPulsedResource_v2(world, p))
      world$resource <- Resource.CC
    }
    
    if(input$resource == "resources_drifting"){
      par0 <- getCCpars(mu_x0 = as.numeric(input$mu.x0), 
                        mu_t0 = as.numeric(input$mu.t0),
                        beta_x = as.numeric(input$beta.x),
                        beta_t = as.numeric(input$beta.t),
                        n.years = as.numeric(input$years),
                        sigma_x = as.numeric(input$x.sd),
                        sigma_t = as.numeric(input$t.sd))
      
      Resource.CC <- aaply(par0, 1, function(p) getPulsedResource(world, p))
      world$resource <- Resource.CC
    }
    
    parameters <- c(epsilon = as.numeric(input$epsilon), 
                    alpha = as.numeric(input$alphabeta) * as.numeric(input$kappa),
                    beta = as.numeric(input$alphabeta) * (1-as.numeric(input$kappa)),
                    gamma = as.numeric(input$gamma))
    
    sim <- runManyYears(World=world, parameters = parameters, 
                        n.years = as.numeric(input$years), 
                        threshold = as.numeric(input$threshold), verbose=FALSE)
    
    
    indices <- data.frame(computeIndices(sim[[length(sim)]], 
                                         world$resource[length(sim)-1,,], world), 
                          final_similarity = computeEfficiency(sim[[length(sim)-1]], 
                                                               sim[[length(sim)]], world))
    indices <- format(indices, nsmall=4)
    newlist <- list(sim,indices)
    
  })
  
  resourceImage <- eventReactive(input$run,{
    if(input$world == "world_optimal"){
      world <- getOptimalPop(tau=100, X.min = 0, X.max = 100, dx=1, 
                             x.peak=as.numeric(input$mu.x0), t.peak=as.numeric(input$mu.t0), 
                             x.sd=as.numeric(input$x.sd), t.sd=as.numeric(input$t.sd))
    }
    else{
      data(world)
      world <- get(input$world)
    } 
    
    if(input$resource == "resources_island"){
      par0 <- getCCpars(mu_x0 = as.numeric(input$mu.x0), 
                        mu_t0 = as.numeric(input$mu.t0),
                        beta_x = as.numeric(input$beta.x),
                        beta_t = as.numeric(input$beta.t),
                        n.years = as.numeric(input$years),
                        sigma_x = as.numeric(input$x.sd),
                        sigma_t = as.numeric(input$t.sd))
      
      Resource.CC <- aaply(par0, 1, function(p) getPulsedResource_v2(world, p))
    }
    
    if(input$resource == "resources_drifting"){
      par0 <- getCCpars(mu_x0 = as.numeric(input$mu.x0), 
                        mu_t0 = as.numeric(input$mu.t0),
                        beta_x = as.numeric(input$beta.x),
                        beta_t = as.numeric(input$beta.t),
                        n.years = as.numeric(input$years),
                        sigma_x = as.numeric(input$x.sd),
                        sigma_t = as.numeric(input$t.sd))
      
      Resource.CC <- aaply(par0, 1, function(p) getPulsedResource(world, p))
    }
    
    par(mfrow = c(ceiling(min(dim(Resource.CC))/5), 5), mar = c(0,0,3,0), oma = c(2,2,4,2), tck = 0.01)
    for (i in 1:min(dim(Resource.CC))) image(1:100, 1:100, Resource.CC[i,,], main = paste("year", i-1), yaxt = "n", xaxt = "n")
    
  })
  
  output$Image <- renderPlot({
    plotManyRuns(simulation()[[1]], nrow=ceiling(length(simulation()[[1]])/6), labelyears=TRUE)
  }, res = 150)
  
  output$Resourceimage <- renderPlot({
    resourceImage()
  }, res = 150)
  
  output$Indices <- renderTable({
    simulation()[[2]]
  })
}
