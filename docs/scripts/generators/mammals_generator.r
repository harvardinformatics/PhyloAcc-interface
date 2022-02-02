############################################################
# For turtles web, 07.21
# This generates the file "mammals.html"
############################################################


cat("Rendering mammals.rmd/html\n")
Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/RStudio/bin/pandoc/")
library(rmarkdown)
setwd("C:/bin/PhyloAcc-interface/docs/scripts/generators/")
output_dir = "../.."
render("../markdown/mammals.rmd", output_dir = output_dir, params = list(output_dir = output_dir), quiet = TRUE)