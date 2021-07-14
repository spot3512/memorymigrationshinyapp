pcks <- c("shiny","memorymigration", "shinythemes", "rsconnect")
lapply(pcks, library, character = TRUE)

library("shiny")
library("memorymigration")
library("shinythemes")
library("rsconnect")

# Setting up Input Interface --------------------------
ui <- fluidPage(
  theme = shinytheme("superhero"),
  h1("Memory Migration model"),
  fluidRow(
    column(3,
           h3("Migratory Population"),
           actionButton(inputId = "run", 
                        label = "Run Model"),
           downloadButton(outputId = "downloadData", 
                          label = "Download Results"),
           #   progressBar(id = "pb1", value = 72),
           radioButtons(inputId = "world",
                        label = "Initial distribution of population in year 0", 
                        choices = c("Optimal" = "world_optimal", "Non-Migratory" = "world_nonmigratory", "Sinusoidal" = "world_sinusoidal"), inline = TRUE),
           numericInput(inputId = "years", label = "Duration of simulation:", value = 5, step = 1),
           numericInput(inputId = "threshold", label = "Threshold of Similarity", value = 0.999, min = 0, max = 1, step = 0.01),
           numericInput(inputId = "alpha",
                        label = "Resource Following - Alpha",
                        value = 100, min = 0, step = 1),
           numericInput(inputId = "beta",
                        label = "Strength of Sociality - Beta",
                        value = 400, min = 0, step = 72),
           sliderInput(inputId = "kappa",
                       label = "Proportion reference vs. working memory - Kappa",
                       value = 0, min = 0, max = 1),
           numericInput(inputId = "lambda",
                        label = "Spatial Scale of Sociality - Lambda",
                        value = 80, min = 0, step = 1),
           numericInput(inputId = "epsilon",
                        label = "Diffusion Parameter - Epsilon",
                        value = 4, min = 0, step = 1)
    ),
    column(6,
           h3("Migratory Population"),
           tableOutput("Indices"),
           plotOutput("Image", height = "400px"),
           h3("Memory"),
           plotOutput("Memory", height = "400px"),
           plotOutput("MigrationHat", height = "400px"),
           h3("Resource"),
           plotOutput("Resourceimage", height = "400px")
    ),
    column(3,
           h3("Resource"),
           actionButton(inputId = "viewresource", 
                        label = "View Resource"),
           sliderInput(inputId = "x.sd",
                       label = "Resource Space Distribution",
                       value = 12, min = 0, max = 15),
           sliderInput(inputId = "t.sd",
                       label = "Resource Time Distribution",
                       value = 6, min = 0, max = 15),
           numericInput(inputId = "mu.x0", label = "Initial Resource Position in Space", value = 30, step = 1),
           numericInput(inputId = "mu.t0", label = "Initial Resource Position in Time", value = 25, step = 1),
           numericInput(inputId = "beta.x", label = "Resource Change in Space", value = 0, step = 1),
           numericInput(inputId = "beta.t", label = "Resource Change in Time ", value = 0, step = 1),
           numericInput(inputId = "psi_x", label = "Stochasticity in Space", value = 0, step = 1),
           numericInput(inputId = "psi_t", label = "Stochasticity in Time ", value = 0, step = 1),
           numericInput(inputId = "x_null", label = "Constraint in Space ", value = 72, step = 1),
           numericInput(inputId = "n_years_null", label = "Years of Stable Resource ", value = 0, step = 1),
           radioButtons(inputId = "resource",
                        label = "Type of resource", 
                        choices = c("Island" = "resources_island", "Drifting" = "resources_drifting"), inline = TRUE)
    ),
    
    
  ))



server <- function(input, output, session) {
  #pcks <- c("shiny","sf","ggplot2","magrittr","plyr", "gplots", 
  #          "memorymigration", "DT", "ggthemes", "minpack.lm", "fields","scales")
  #lapply(pcks, require, character = TRUE)
  
  
  
  # Setting up World ---------------------
  
  simulation <- eventReactive(input$run, {
    
    if(input$world == "world_optimal"){
      world <- getOptimalPop(tau=100, X.min = -100, X.max = 100, dx=1, 
                             x1 =as.numeric(input$mu.x0), 
                             x2 = -as.numeric(input$mu.x0),
                             t.peak=as.numeric(input$mu.t0), 
                             x.sd=as.numeric(input$x.sd), 
                             t.sd=as.numeric(input$t.sd))
    }
    
    if(input$world == "world_nonmigratory"){
      world <- getSinePop(tau = 100, peak.max = 1, peak.min = -1, sd = 10)
    }
    if(input$world == "world_sinusoidal"){
      world <- getSinePop(tau = 100, peak.max = as.numeric(input$mu.x0), peak.min = -as.numeric(input$mu.x0), sd = 10)
    }
    world$m0 <- fitMigration(t = world$time, x = getMem(world$pop, world))
    
    par0 <- getCCpars(mu_x0 = as.numeric(input$mu.x0), 
                      mu_t0 = as.numeric(input$mu.t0),
                      beta_x = as.numeric(input$beta.x),
                      beta_t = as.numeric(input$beta.t),
                      n.years = as.numeric(input$years),
                      sigma_x = as.numeric(input$x.sd),
                      sigma_t = as.numeric(input$t.sd),
                      psi_x = as.numeric(input$psi_x), 
                      psi_t = as.numeric(input$psi_t),
                      n.years.null = as.numeric(input$n_years_null))
    if(input$resource == "resources_island"){ 
      Resource.CC <- aaply(par0, 1, function(p) getResource_island(world, p))
      world$resource <- Resource.CC
    }
    
    if(input$resource == "resources_drifting"){
      Resource.CC <- aaply(par0, 1, function(p) getResource_drifting(world, p, as.numeric(input$x_null)))
      world$resource <- Resource.CC
    }
    attr(world$resource, "par") <- par0[nrow(par0),]
    resource_param <- data.frame(mu_x0 = as.numeric(input$mu.x0), 
                                 mu_t0 = as.numeric(input$mu.t0),
                                 beta_x = as.numeric(input$beta.x),
                                 beta_t = as.numeric(input$beta.t),
                                 n.years = as.numeric(input$years),
                                 sigma_x = as.numeric(input$x.sd),
                                 sigma_t = as.numeric(input$t.sd),
                                 psi_x = as.numeric(input$psi_x), 
                                 psi_t = as.numeric(input$psi_t),
                                 n.years.null = as.numeric(input$n_years_null),
                                 world = input$world,
                                 resource = input$resource) 
    
    parameters <- c(epsilon = as.numeric(input$epsilon), 
                    alpha = as.numeric(input$alpha),
                    beta = as.numeric(input$beta),
                    kappa = as.numeric(input$kappa),
                    lambda = as.numeric(input$lambda))
    
    param.df <- data.frame(
      epsilon = as.numeric(input$epsilon), 
      alpha = as.numeric(input$alpha),
      beta = as.numeric(input$beta),
      kappa = as.numeric(input$kappa),
      lambda = as.numeric(input$lambda)
    )
    
    ## Running the model ---------------------
    sim <- runManyYears(world=world, parameters = parameters, 
                        n.years = as.numeric(input$years), 
                        threshold = as.numeric(input$threshold), 
                        verbose=TRUE)
    
    indices <- data.frame(computeIndices(sim$pop[[length(sim)]], 
                                         world$resource[length(sim$pop)-1,,], world),
                          avgFE = computeAvgEfficiency(sim$pop, world$resource, world),
                          TE = computeTotalError(sim, world),
                          avgTE = computeAvgTotalError(sim, world, par0),
                          final_similarity = computeEfficiency(sim$pop[[length(sim$pop)-1]], 
                                                               sim$pop[[length(sim$pop)]], world), 
                          n.runs = length(sim$pop) - 1,
                          resource_param, param.df)
    if(as.numeric(input$beta.x) != 0){ 
      indices$SA_total <- computeSpatialAdaptationIndex(sim, resource_param)
    }else {indices$SA_total <- NA}
    
    
    #parameters.df <- ldply (parameters, data.frame)
    indices <- format(indices, digits=4)
    
    ## Plotting migrations -----------------
    
    list(sim = sim,
         indices = indices,
         world = world, 
         x.peak = as.numeric(input$x.peak),
         t.peak = as.numeric(input$t.peak))
  })
  
  ## Resource Image ---------------------
  resourceImage <- eventReactive(input$run | input$viewresource,{
    if(input$world == "world_optimal"){
      world <- getOptimalPop(tau=100, X.min = -100, X.max = 100, dx=1, 
                             x1 =as.numeric(input$mu.x0), 
                             x2 = -as.numeric(input$mu.x0),
                             t.peak=as.numeric(input$mu.t0), 
                             x.sd=as.numeric(input$x.sd), 
                             t.sd=as.numeric(input$t.sd))
    }
    if(input$world == "world_nonmigratory"){
      world <- getSinePop(tau = 100, peak.max = 1, peak.min = -1, sd = 10)
    }
    if(input$world == "world_sinusoidal"){
      world <- getSinePop(tau = 100, peak.max = as.numeric(input$mu.x0), peak.min = -as.numeric(input$mu.x0), sd = 10)
    }
    world$m0 <- fitMigration(t = world$time, x = getMem(world$pop, world))
    par0 <- getCCpars(mu_x0 = as.numeric(input$mu.x0), 
                      mu_t0 = as.numeric(input$mu.t0),
                      beta_x = as.numeric(input$beta.x),
                      beta_t = as.numeric(input$beta.t),
                      n.years = as.numeric(input$years),
                      sigma_x = as.numeric(input$x.sd),
                      sigma_t = as.numeric(input$t.sd),
                      psi_x = as.numeric(input$psi_x), 
                      psi_t = as.numeric(input$psi_t),
                      n.years.null = as.numeric(input$n_years_null))
    if(input$resource == "resources_island"){ 
      Resource.CC <- aaply(par0, 1, function(p) getResource_island(world, p))
      world$resource <- Resource.CC
    }
    
    if(input$resource == "resources_drifting"){
      Resource.CC <- aaply(par0, 1, function(p) getResource_drifting(world, p, as.numeric(input$x_null)))
      world$resource <- Resource.CC
    }
    attr(world$resource, "par") <- par0[nrow(par0),]
    
    
    
    par(mfrow = c(ceiling(min(dim(Resource.CC))/5), 5), mar = c(1,1,1,1), oma = c(2,2,0,2), tck = 0.01)
    for (i in 1:min(dim(Resource.CC))) image(Resource.CC[i,,], main = paste("year", i-1), yaxt = "n", xaxt = "n")
    
  })
  
  output$Image <- renderImage({
    
    
    progress <- Progress$new(session, min=1, max=15)
    on.exit(progress$close())
    
    progress$set(message = 'Calculation in progress',
                 detail = 'This may take a while...')
    
    width  <- session$clientData$output_Image_width
    height <- session$clientData$output_Image_height
    
    pixelratio <- session$clientData$pixelratio
    
    outfile <- tempfile(fileext='.png')
    
    png(outfile, width = width*pixelratio, height = height*pixelratio,
        res = 72*pixelratio)
    plotManyRuns(simulation()[[1]]$pop, world = simulation()[[3]], nrow=ceiling(length(simulation()[[1]]$pop)/10), labelyears=TRUE)
    dev.off()
    
    list(src = outfile,
         width = width,
         height = height,
         alt = "This is alternate text")
  }, deleteFile = TRUE)
  
  
  output$Resourceimage <- renderImage({
    width  <- session$clientData$output_Resourceimage_width
    height <- session$clientData$output_Resourceimage_height
    
    pixelratio <- session$clientData$pixelratio
    
    outfile <- tempfile(fileext='.png')
    
    png(outfile, width = width*pixelratio, height = height*pixelratio,
        res = 72*pixelratio)
    resourceImage()
    dev.off()
    
    list(src = outfile,
         width = width,
         height = height,
         alt = "This is alternate text")
  }, deleteFile = TRUE)
  
  output$Indices <- renderTable({
    simulation()[[2]][c("FE", "avgFE", "TE", "avgTE", "SA_total", "final_similarity", "n.runs" )]
    
  }, digits = 3)
  
  
  # double Plot (migration and resource) ----------------------
  output$Memory <- renderImage({
    width  <- session$clientData$output_Memory_width
    height <- session$clientData$output_Memory_height
    
    pixelratio <- session$clientData$pixelratio
    
    outfile <- tempfile(fileext='.png')
    
    png(outfile, width = width*pixelratio, height = height*pixelratio,
        res = 72*pixelratio)
    doublePlotForShiny(simulation()$sim$pop, 
                       world = simulation()$world)
    dev.off()
    
    list(src = outfile,
         width = width,
         height = height,
         alt = "This is alternate text")
  }, deleteFile = TRUE)
  
  output$MigrationHat <- renderImage({
    width  <- session$clientData$output_MigrationHat_width
    height <- session$clientData$output_MigrationHat_height
    
    pixelratio <- session$clientData$pixelratio
    
    outfile <- tempfile(fileext='.png')
    
    png(outfile, width = width*pixelratio, height = height*pixelratio,
        res = 72*pixelratio)
    par(mfrow = c(1,2), mar = c(2,2,2,2), 
        tck = 0.01, mgp = c(1.5,.25,0), 
        bty = "l", cex.lab = 1.25, las = 1, xpd = NA)
    with(simulation(),
         plotMigrationHat(sim$migration.hat, 
                          x.peak = x.peak, t.peak = t.peak)
    )
    dev.off()
    
    list(src = outfile,
         width = width,
         height = height,
         alt = "This is alternate text")
  }, deleteFile = TRUE)
  
  output$downloadData <- downloadHandler(
    filename = "simulationRun.csv",
    content = function(file) {
      
      write.csv(simulation()[[2]], file)
    }
  )
}



shinyApp(ui, server)

