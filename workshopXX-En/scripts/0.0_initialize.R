# Description  ------------------------------------------------------------
##########################################################################################
# Created by Marc-Olivier Beausoleil
# McGill University 
# Created Tuesday, February 12, 2020 
# Why: 
# Requires 
# NOTES: 
  # loading a bunch of packages efficiently https://statsandr.com/blog/an-efficient-way-to-install-and-load-r-packages/
##########################################################################################

# Load libraries ----------------------------------------------------------
librarian::shelf(tidyverse, plotly, viridis, purrr, dplyr, 
                 ggplot2, ggExtra, gridExtra, 
                 gavinsimpson / ggvegan, 
                 scales, MASS, Rlab, fitdistrplus, dummies, 
                 pwr, lme4, broom.mixed, lmerTest, 
                 sf, mapview, 
                 cds, TeachingDemos, remotes)

# Load packages and install  ----------------------------------------------
# # Package names
# packages.all <- c("tidyverse", "plotly", "viridis", "purrr", "dplyr", 
#                   "ggplot2", "ggExtra", "gridExtra", 
#                   # "ggvegan", 
#                   "gavinsimpson / ggvegan", 
#                   "scales", "MASS", "Rlab", "fitdistrplus", "dummies", 
#                   "pwr", "lme4", "broom.mixed", "lmerTest", 
#                   "sf", "mapview", 
#                   "cds", "TeachingDemos", "remotes")
# 
# # Check which package is installed 
# installed_packages <- packages.all %in% rownames(installed.packages())
# 
# # Install the pacakges on GitHub first
# g.veg.pos = which(packages.all == "ggvegan")
# if (grepl(pattern = "ggvegan",packages.all[!installed_packages])) {
#   # ggvegan is only on github: 
#   remotes::install_github("gavinsimpson/ggvegan")
# }
# 
# # Clear GitHub package from list 
# packages = packages.all[-g.veg.pos]
# 
# # Get packages on Cran 
# installed_packages <- packages.all %in% rownames(installed.packages())
# 
# # Install on Cran 
# if (any(installed_packages == FALSE)) {
#   install.packages(packages[!installed_packages])
# }
# 
# # Load packages
# invisible(lapply(packages, library, character.only = TRUE))


# Using pacman ------------------------------------------------------------
# pacman::p_load(tidyverse,plotly,viridis,purrr,dplyr,
#                ggplot2,ggExtra,gridExtra,ggvegan,
#                scales,MASS,Rlab,fitdistrplus,dummies,
#                pwr,lme4,broom.mixed,lmerTest,
#                sf,mapview,
#                cds,TeachingDemos, remotes)

# old school --------------------------------------------------------------
# library(tidyverse)
# library(plotly)
# library(viridis)
# library(purrr) 
# library(dplyr)
# library(ggplot2)
# library(ggExtra)
# library(gridExtra)
# library(ggvegan)
# library(scales)
# library(MASS) # mvrnorm
# library(Rlab) # For distribution Rlab::dbern()
# library(fitdistrplus) # descdist and fitdist()
# library(dummies)
# library(pwr)
# library(lme4,warn.conflicts = FALSE)
# library(broom.mixed)
# library(lmerTest, warn.conflicts = F)
# library(sf, warn.conflicts = FALSE)
# library(mapview, warn.conflicts = FALSE)
# library(cds)
# library(TeachingDemos)


# suppressPackageStartupMessages( library(dplyr) )

# Source scripts ----------------------------------------------------------
source(file = "scripts/Galton_board.R")
source(file = "scripts/marginal_plot.R")


# Get session information -------------------------------------------------
dir.create("workshopXX-En/output/session_info",recursive = TRUE, showWarnings = FALSE)
sink("workshopXX-En/output/session_info/session_information.txt",append = FALSE)
# sink("~/Desktop/session_information.txt",append = FALSE)
cat("##### R Version Information ############################################################\n")
version

cat("\n\n##### Collect Information About the Current R Session ############################################################\n\n")
sessionInfo() 
# loadedNamespaces()
sink()

cat("Done!\n")


