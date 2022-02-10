# Description  ------------------------------------------------------------
#### ### ### ## #### ### ### ## #### ### ### ## 
# Extract the script from the RMarkdown 
# Created by Marc-Olivier Beausoleil 
# 
# Why: 
# Requires:
# NOTES: 
# Reference : 
#### ### ### ## #### ### ### ## #### ### ### ## 

library(knitr)
file.exists("workshopXX-En/workshopXX-en.Rmd")
knitr::purl(input = "workshopXX-En/workshopXX-en.Rmd", output = "workshopXX-En/scripts/workshopXX-en.r")
