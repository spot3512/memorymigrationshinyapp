pcks <- c("shiny","sf","ggplot2","magrittr","plyr", "gplots", "memorymigration")
lapply(pcks, require, character = TRUE)

data(world)
data(resources_island)
data(resources_drifting)

source("ShinyScripts.R")
shinyApp(ui, server)
