# Description  ------------------------------------------------------------
##########################################################################################
# Created by Marc-Olivier Beausoleil
# McGill University 
# Created Tuesday, February 12, 2020 
# Why: 
# Requires 
  if("librarian" %in% rownames(installed.packages()) == FALSE) {install.packages("librarian")}
  # "scripts/Galton_board.R"
  # "scripts/marginal_plot.R"
  # Install "ImageMagick" https://imagemagick.org https://imagemagick.org/script/download.php
# NOTES: 
  # loading a bunch of packages efficiently https://statsandr.com/blog/an-efficient-way-to-install-and-load-r-packages/
##########################################################################################

# Load libraries ----------------------------------------------------------
librarian::shelf(tidyverse,
                 broom.mixed,
                 dplyr, DT,
                 cds, 
                 # dummies, # not availabie anymore https://cran.r-project.org/web/packages/dummies/index.html
                 fontawesome,
                 fitdistrplus, 
                 gavinsimpson / ggvegan, 
                 ggplot2, ggExtra, gridExtra, 
                 lme4, lmerTest, 
                 mapview, MASS, 
                 plotly, purrr, pwr, 
                 remotes, Rlab, 
                 scales, shiny, sf, 
                 TeachingDemos, 
                 viridis, 
                 cran_repo = 'https://cran.r-project.org')
# suppressPackageStartupMessages( library(dplyr) )

# Source scripts ----------------------------------------------------------
source(file = "scripts/Galton_board.R")
source(file = "scripts/marginal_plot.R")
source(file = "scripts/draw.normal.R")

# Get session information -------------------------------------------------
dir.create("output/session_info",recursive = TRUE, showWarnings = FALSE)
sink("output/session_info/session_information.txt",append = FALSE)
# sink("~/Desktop/session_information.txt",append = FALSE)
cat("##### R Version Information ############################################################\n")
version

cat("\n\n##### Collect Information About the Current R Session ############################################################\n\n")
sessionInfo() 
# loadedNamespaces()
sink()

# cat("Done!\n")
