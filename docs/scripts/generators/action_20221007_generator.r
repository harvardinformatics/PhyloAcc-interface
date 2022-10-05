############################################################
# For turtles web, 07.21
# This generates the file "action_20221007.html"
############################################################


cat("Rendering action_20221007.rmd/html\n")

library(here)
setwd(here("docs", "scripts", "generators"))

Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/RStudio/bin/pandoc/")
library(rmarkdown)

output_dir = "../.."
render("../markdown/action_20221007.rmd", output_dir = output_dir, params = list(output_dir = output_dir), quiet = TRUE)