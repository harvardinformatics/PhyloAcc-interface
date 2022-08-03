############################################################
# For turtles web, 07.21
# This generates the file "sim_trees.html"
############################################################


cat("Rendering sim_trees.rmd/html\n")

library(here)
setwd(here("docs", "scripts", "generators"))

Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/RStudio/bin/pandoc/")
library(rmarkdown)

output_dir = "../.."
render("../markdown/sim_trees.rmd", output_dir = output_dir, params = list(output_dir = output_dir), quiet = TRUE)