# Description  ------------------------------------------------------------
##########################################################################################
# Created by Marc-Olivier Beausoleil
# McGill University 
# Created ~10 June 2022
# Why: 
# Requires 
# NOTES: 
# devtools::install_github("jsugarelli/packagefinder")
# https://github.com/jsugarelli/packagefinder
##########################################################################################

library(packagefinder)
# if return.df = TRUE, you can play with the df.
sim.pack = findPackage(c("simulation","simulate", "DNA"), mode = "or", display = "console",return.df = FALSE, limit.results = 100)
sim.pack = findPackage(c("site occupancy"), mode = "and", display = "console",return.df = FALSE, limit.results = 100)
sim.pack = findPackage(c("simulation", "DNA"), mode = "and", display = "console",return.df = FALSE, limit.results = 100)
sim.pack = findPackage(c("simulate", "DNA"), mode = "and", display = "console",return.df = FALSE, limit.results = 100)
findPackage(c("simulation","simulate"), mode = "or", display = "console")
# findPackage(c("simulation","simulate"), mode = "or", display = "browser",return.df = FALSE)
nrow(sim.pack)
names(sim.pack)
View(sim.pack[,c("Name","Short Description")])

install.packages("fwsim")
library(fwsim)
example(fwsim)

install.packages("tidydice")
library(tidydice)
roll_dice(times = 10, rounds = 10) %>% 
  plot_dice()

roll_dice(sides = 4,times = 5, rounds = 2) %>% 
  plot_dice()
flip_coin(times = 2,rounds = 5) %>% 
  plot_dice()
sample(0:1, 2, replace = T)