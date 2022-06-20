## ----setup, echo=F, warning=FALSE-----------------------------------------------------------------------------------------------------------------------------
# in the end, should be about 1500 lines
knitr::opts_chunk$set(
  comment = "#",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  fig.width=6, fig.height=6,
  fig.retina = 3,
  fig.align = 'center'
)

if("librarian" %in% rownames(installed.packages()) == FALSE) {install.packages("librarian")}
library(librarian) # even better than pacman 
# library(pacman)
# Need to install all packages in initialize
suppressWarnings(source(file = "scripts/0.0_initialize.R"))
# Check also suppressMessages()


## ----normal_compare_theoretical_simulated, echo=FALSE, fig.width=8, fig.height=4------------------------------------------------------------------------------
par(mfrow=c(1,2), mar = c(3,4,.5,.5), cex = 1.3)
set.seed(12345)
x = seq(-5,5,by=.1)
y = dnorm(x)

curve(expr = dnorm(x), 
      from = -5,
      to = 5, 
      ylim = c(0,1), 
      xlab = "", 
      ylab = "Density", 
      lwd =3);abline(v = 0, lty =3)
  title(xlab="x", line=2, cex.lab=1.2)
polygon(c(x, rev(x), 0), c(y, rep(0,length(y)), 0), col=scales::alpha("blue",.5))
x2 = rnorm(6)
y2 = density(x2)
plot(y2, 
     main = "", 
     xlab = "", 
     ylab = "Density", 
     xlim = c(-5,5), 
     ylim = c(0,1),
     lwd = 3);abline(v = 0, lty =3)
  title(xlab="x", line=2, cex.lab=1.2)
polygon(c(y2$x, rev(y2$x), 0), c(y2$y, rep(0,length(y2$y)), 0), col=scales::alpha("blue",.5))



## ----r_tips2, echo=FALSE, fig.width=13,fig.height=6.5---------------------------------------------------------------------------------------------------------
# Description  ------------------------------------------------------------
#### ### ### ## #### ### ### ## #### ### ### ## 
# Genetic drift simulation
# Created by Marc-Olivier Beausoleil
# 2022-01-07
# Why: 
# Requires:
# NOTES: 
# Drift is (from Futuyma)
# - unbiased
# - random fluctuations in allele frequency are larger in smaller populations
# - drift causes genetic variation to be lost
# - drift causes populations that are initially identical to become different
# - an allele can become fixed without the benefit of natural selection
# Reference : 
# Futuyma p. 167, figure 7.2
#### ### ### ## #### ### ### ## #### ### ### ## 

# graphing parameters -----------------------------------------------------
par(mfrow = c(2,2), mar = c(3,3,2,2), cex = 1.2)
# Random seed -------------------------------------------------------------
set.seed(1245)

# Simulation parameters ---------------------------------------------------
# Number of gametes to chose from 
n.sperm = 2
n.eggs = 2
# Number of generations (x axis)
gen = 500
# Number of replicate populations 
popu = 5

# variance = p*(1-p)/(2*N)
# Variation is smaller when the population size is bigger 

# Loops -------------------------------------------------------------------
# Loops for all population replicates, tracking allele frequency change over the generations 
# number of individual per population 
n.id.pop = 5*10^seq(0,3, by=1)
# Loop that will change the maximum number of individual per population 
for (l in n.id.pop) {
  # Initial allele frequency 
  p.init = .5 
  # Maximum population size 
  max.pop = l
  # Total number of gametes in the population 
  n.gametes = c(max.pop*(n.sperm+n.eggs))
  # Make an empty object to record the population information 
  all.pops = NULL
  # Loop to track the all population allele frequency change 
  for (j in 1:popu) {
    all.fq.change = .5
    # Loop to track the within population allele frequency change 
    for (i in 1:gen) {
      # If the first iteration, make the probability equal the initial allele frequency 
      if (i == 1) {prob.p = p.init} else {prob.p = prop.all[2]}
      # binomial function to generate the new allele frequency (0 = q, 1 = p)
      allele.fq = rbinom(n = n.gametes, size = 1, prob = prob.p)
      # Randomly sample the population (this is the drift, a random sample of the population)
      all.drift = sample(x = allele.fq, size = max.pop, replace = F)
      # Get the proportion of the alleles in the new population 
      prop.all = prop.table(table(all.drift))
      # Record the p allele only 
      all.fq.change = c(all.fq.change, prop.all[2])
      # If there is an allele that goes to fixation, it'll print NA. In this case, break the for loop and go to the next iteration
      if(is.na(prop.all[2])) {break}
    } # End i
    # Record all population information
    one.pop = data.frame(p.fq = as.numeric(all.fq.change), pop = j)
    all.pops = rbind(all.pops,one.pop)
    # If an allele is fixated, it's going to be recorded 
    all.pops[is.na(all.pops)] <- ifelse(names(prop.all)=="1",yes = 1,no = 0)

  } # End j
  
  # Plot --------------------------------------------------------------------
  # Make the empty plot 
  plot(all.pops$p.fq~c(1:nrow(all.pops)), 
       col = as.factor(all.pops$pop),
       main = paste0("Pop N=",max.pop, ", Start p=",p.init),
       ylab = "",
       xlab = "",
       ylim = c(0,1),
       xlim = c(1,gen),
       type = "n")
  if (l%in% n.id.pop[c(1,3)]) {
  title(ylab="Allele fq p", line=2, cex.lab=1.2)
  }
  if (l%in% n.id.pop[c(3,4)]) {
  title(xlab="Generations", line=2, cex.lab=1.2)
  }
  
  
  # Add the lines per population and colour them 
  for (k in 1:popu) {
    pttmp = all.pops[all.pops$pop==k,]
    points(c(1:nrow(pttmp)), pttmp$p.fq, 
           type = "l",
           col = k)
    final.allele.fq = ifelse(pttmp$p.fq[length(pttmp$p.fq)] %in% c(0),"red",
                             ifelse(pttmp$p.fq[length(pttmp$p.fq)] %in% c(1),"green","black"))
    points(c(nrow(pttmp)), pttmp$p.fq[length(pttmp$p.fq)], 
           type = "p",
           col = scales::alpha(final.allele.fq, .5), pch = 19)
           # col = scales::alpha(c("black","red","green","blue","cyan")[k], .5), pch = 19)
  } # End k
} # End l



## ----fake_fitness_functions, echo=FALSE, fig.width=8,fig.height=3---------------------------------------------------------------------------------------------
# Fake fitness landscapes 
# Fitness test theoretical fitness landscapes 
par(mfrow=c(1,3))
q.fun = function(x, 
                 exponent = 2,
                 factor = -1,
                 xlab = xlab,
                 ylab = "Fitness",
                 lwd = 3,
                 ylim = c(0,4000),
                 cex.text = 1,
                 col.main = "black", 
                 col.line = "red", 
                 col.box = "black", 
                 col.text = "black") {
  y = x^exponent*factor
  plot(y~x, type = "l", 
       axes=FALSE,
       xaxs="i",yaxs="i",
       frame.plot=FALSE, 
       xlab="", ylab = "",lwd = lwd, col = col.main,
       # ylim = c(min(y),max(y)+5000)
       ylim = ylim)
  box(bty="l", lwd = 3)
  mtext(side=2,text=ylab,line = 1.1, col = col.text, cex = cex.text)
  mtext(side=1,text=xlab,line = 1.6, col = col.text, cex = cex.text)
  # axis(side = 1, labels = FALSE, tck = 0.000000001, lwd = 4, col = col.box)
  # axis(side = 2, labels = FALSE, tck = 0.000000001, lwd = 4, col = col.box)
  # y2 = -3*x^2+2500
  # points(y2~x, type = "l", lty = 2, lwd = lwd,
  #        col = col.line)
}
# col.main = "white"
col.main = "black"
xlab = ""
q.fun(-100:100,exponent = 1,factor = 10,xlab = "(Directional)",lwd = 5,cex.text = 2,
      ylab = "Fitness", col.line = NA,ylim = c(-1000, 2000),
      col.main = col.main,col.box = col.main, col.text = col.main)
q.fun(-100:100,xlab = "(Stabilizing)",lwd = 5,cex.text = 2,
      ylab = "", col.line = NA,ylim = c(-5000, 2000),
      col.main = col.main,col.box = col.main, col.text = col.main)
q.fun(-70:70,factor = 1,xlab = "(Disruptive)",lwd = 5,cex.text = 2,
      ylab = "", col.line = NA,ylim = c(-200, 7000),
      col.main = col.main,col.box = col.main, col.text = col.main)



## ----normalX_Y, echo=FALSE, fig.show='hide', fig.width=8,fig.height=4-----------------------------------------------------------------------------------------
# source(file = "scripts/marginal_plot.R")
set.seed(123)
# par(mfrow = c(1,2), mar =c(4,4,3,3), cex = 1.4)
n = 250
x.1 = rnorm(n, mean = 15, sd = 5)
y.1 = 2*x.1 +rnorm(n, mean = 0, sd = 4)#rnorm(n, mean = 5, sd = 2)

x.2 = rnorm(n, mean = 15, sd = 1)
y.2 = 2*x.2 +rnorm(n, mean = 0, sd = .5) # rnorm(n, mean = 5, sd = .5)
# marginal_plot(x.1,y.1,ylim = range(c(y.1,y.2)), xlim = range(c(x.1,x.2)), pch = 19, col = scales::alpha("black",.8), ylab = "Y",xlab = "X")
# marginal_plot(x.2,y.2,ylim = range(c(y.1,y.2)), xlim = range(c(x.1,x.2)), pch = 19, col = scales::alpha("black",.8), ylab = "Y",xlab = "X")

library(ggplot2)
# Use base R here 
# https://stackoverflow.com/questions/71052975/how-to-plot-histograms-in-base-r-in-the-margin-of-a-plot?noredirect=1#comment125604177_71052975 
library(ggExtra)

size.line = .8
text.size = 16
col.line = alpha("black",.5)
df.1 <- data.frame(x = x.1, y = y.1)
df.2 <- data.frame(x = x.2, y = y.2)
p.1 <- ggplot(df.1, aes(x, y)) + 
  geom_point(colour = alpha("black",.5)) + 
  lims(x = range(c(x.1,x.2)), y = range(c(y.1,y.2))) +
  theme_classic() + 
  theme(axis.ticks = element_line(colour = "black"),
        axis.title = element_text(size = text.size),
        axis.text = element_text(size = text.size, colour = "black"),
        axis.text.x = element_text(size = text.size),
        axis.text.y = element_text(size = text.size))

p.2 <- ggplot(df.2, aes(x, y)) + 
  geom_point(colour = alpha("black",.5)) + 
  lims(x = range(c(x.1,x.2)), y = range(c(y.1,y.2))) +
  theme_classic() + 
  theme(axis.ticks = element_line(colour = "black"),
        axis.title = element_text(size = text.size),
        axis.text = element_text(size = text.size, colour = "black"),
        axis.text.x = element_text(size = text.size),
        axis.text.y = element_text(size = text.size))
p.1
p.2

p.h1 <- p.1 + geom_hline(yintercept=mean(y.1), linetype="dashed", color = col.line, size=size.line)
p.h1
p.h.v1 <- p.h1 + geom_vline(xintercept=mean(x.1), linetype="dashed", color = col.line, size=size.line)
p.h.v1
p.h2 <- p.2 + geom_hline(yintercept=mean(y.2), linetype="dashed", color = col.line, size=size.line)
p.h2
p.h.v2 <- p.h2 + geom_vline(xintercept=mean(x.2), linetype="dashed", color = col.line, size=size.line)
p.h.v2

positions.1 <- data.frame(
  x = c(mean(x.1)-sd(x.1), mean(x.1)+sd(x.1), mean(x.1)+sd(x.1), mean(x.1)-sd(x.1)),
  y = c(mean(y.1)-sd(y.1), mean(y.1)-sd(y.1), mean(y.1)+sd(y.1), mean(y.1)+sd(y.1))
)
positions.2 <- data.frame(
  x = c(mean(x.2)-sd(x.2), mean(x.2)+sd(x.2), mean(x.2)+sd(x.2), mean(x.2)-sd(x.2)),
  y = c(mean(y.2)-sd(y.2), mean(y.2)-sd(y.2), mean(y.2)+sd(y.2), mean(y.2)+sd(y.2))
)

positions.x1 <- data.frame(
  x = c(mean(x.1)-sd(x.1), mean(x.1)+sd(x.1), mean(x.1)+sd(x.1), mean(x.1)-sd(x.1)),
  y = c(-Inf, -Inf, Inf, Inf)
)
positions.x2 <- data.frame(
  x = c(mean(x.2)-sd(x.2), mean(x.2)+sd(x.2), mean(x.2)+sd(x.2), mean(x.2)-sd(x.2)),
  y = c(-Inf, -Inf, Inf, Inf)
)
positions.y1 <- data.frame(
  x = c(Inf, -Inf, -Inf, Inf),
  y = c(mean(y.1)-sd(y.1), mean(y.1)-sd(y.1), mean(y.1)+sd(y.1), mean(y.1)+sd(y.1))
)
positions.y2 <- data.frame(
  x = c(Inf, -Inf, -Inf, Inf),
  y = c(mean(y.2)-sd(y.2), mean(y.2)-sd(y.2), mean(y.2)+sd(y.2), mean(y.2)+sd(y.2))
)
# p.h.v + geom_polygon(data = positions, aes(x = x, y = y),fill = alpha("gray70",.6))
p.h.v.sd.x1 = p.h.v1 + geom_polygon(data = positions.x1, aes(x = x, y = y),fill = alpha("gray70",.4))
p.h.v.sd.x1
p.h.v.sd.x2 = p.h.v2 + geom_polygon(data = positions.x2, aes(x = x, y = y),fill = alpha("gray70",.4))
p.h.v.sd.x2

p.h.v.sd.x.sd.y1 = p.h.v.sd.x1 + geom_polygon(data = positions.y1, aes(x = x, y = y),fill = alpha("gray70",.4))
p.h.v.sd.x.sd.y1
p.h.v.sd.x.sd.y2 = p.h.v.sd.x2 + geom_polygon(data = positions.y2, aes(x = x, y = y),fill = alpha("gray70",.4))
p.h.v.sd.x.sd.y2

p.h.v.sd.x.sd.y.l1 = p.h.v.sd.x.sd.y1 +  geom_smooth(method = "lm", se = TRUE, col =alpha("red",.5))
p.h.v.sd.x.sd.y.l1
p.h.v.sd.x.sd.y.l2 = p.h.v.sd.x.sd.y2 +  geom_smooth(method = "lm", se = TRUE, col =alpha("red",.5))
p.h.v.sd.x.sd.y.l2

p.marg1 = ggExtra::ggMarginal(p.h.v.sd.x.sd.y.l1, type = "histogram", fill = "gray80", col = "gray70")
p.marg1
p.marg2 = ggExtra::ggMarginal(p.h.v.sd.x.sd.y.l2, type = "histogram", fill = "gray80", col = "gray70")
p.marg2
# plot(y.1~x.1, ylim = range(c(y.1,y.2)), xlim = range(c(x.1,x.2)), pch = 19, col = scales::alpha("black",.8), ylab = "Y",xlab = "X")
# plot(y.2~x.2, ylim = range(c(y.1,y.2)), xlim = range(c(x.1,x.2)), pch = 19, col = scales::alpha("black",.8), ylab = "Y",xlab = "X")



## ----normalX_Y1, echo=FALSE, fig.width=5,fig.height=4---------------------------------------------------------------------------------------------------------
p.1


## ----normalX_Y2, echo=FALSE, fig.width=5,fig.height=4---------------------------------------------------------------------------------------------------------
p.2


## ----normalX_Y3, echo=FALSE, fig.width=5,fig.height=4---------------------------------------------------------------------------------------------------------
# p.h1
p.h.v1


## ----normalX_Y4, echo=FALSE, fig.width=5,fig.height=4---------------------------------------------------------------------------------------------------------
# p.h2
p.h.v2



## ----normalX_Y5, echo=FALSE, fig.width=5,fig.height=4---------------------------------------------------------------------------------------------------------
# p.h.v.sd.x1
p.h.v.sd.x.sd.y1


## ----normalX_Y6, echo=FALSE, fig.width=5,fig.height=4---------------------------------------------------------------------------------------------------------
# p.h.v.sd.x2
p.h.v.sd.x.sd.y2


## ----normalX_Y7, echo=FALSE, fig.width=5,fig.height=4---------------------------------------------------------------------------------------------------------
p.h.v.sd.x.sd.y.l1


## ----normalX_Y8, echo=FALSE, fig.width=5,fig.height=4---------------------------------------------------------------------------------------------------------
p.h.v.sd.x.sd.y.l2


## ----normalX_Y9, echo=FALSE, fig.width=5,fig.height=4---------------------------------------------------------------------------------------------------------
p.marg1


## ----normalX_Y10, echo=FALSE, fig.width=5,fig.height=4--------------------------------------------------------------------------------------------------------
p.marg2


## ----poissonX_Y, echo=FALSE, fig.show='hide', fig.width=8,fig.height=4----------------------------------------------------------------------------------------
set.seed(1217)
n = 250
b0 = 0.8
b1 = 0.5
b2 = 0.5
# x1 = rnorm(n, mean = 5, sd = 1) # 
# x1 = runif(n, min = 1, max = 3)
x1 = runif(n, min = 0, max = 6)
# Linear on the log scale 
# https://aosmith.rbind.io/2018/07/18/simulate-poisson-edition/ 
lambda1 = exp(b0 + b1*x1) # This defines the mean at each point 
# Since lambda is the number of events/time*time = nb of event, 
# each x point will change the number of events (y) by drawing a value from a poisson distribution that has a different mean each time
# The MEAN and VARIANCE of the poisson process does NOT need to be an integer. Only the OUPUT of the poisson distribution NEEDS to be an integer.
# Lambda can be changed in multiple ways: the number of events for a period of time OR change the time interval
# If you take the mean and variance of the y1, you WON'T see the so called equivalence between mean and variance of poisson distribution. Because this is the Poisson PROCESS and not the DISTRIBUTION that you are observing
# https://www.probabilitycourse.com/chapter10/10_1_0_basic_concepts.php 
# Keep in mind that a RANDOM PROCESS is a collection of RANDOM VARIABLES "collected".
# In other words, you might measure a random variable that is the outcome of a random process that you try to model https://web.ma.utexas.edu/users/mks/M358KInstr/RandomVariables.pdf 

# mean(lambda1)
y1 = rpois(n, lambda = lambda1) 
# plot(density((y1)))
# hist(((y1)))
# mean((y1))

plot(log(y1)~x1)
plot(y1~x1);abline(h = seq(range(y1)[1],range(y1)[2],by = 1),lty = 3)
x2 = rnorm(n, mean = 3, sd = 1)

# x2 = runif(n, min = 2, max = 5)
lambda2 = exp(b0 + b2*x2)
# mean(lambda2)
y2 = rpois(n, lambda = lambda2) 
# mean((y2))

library(ggplot2)
library(ggExtra)
library(gridExtra)

size.line = .8
text.size = 16
col.line = alpha("black",.5)
df.1 <- data.frame(x = x1, y = y1, ly = log(y1), mod = "model1")
df.2 <- data.frame(x = x2, y = y2, ly = log(y2), mod = "model2")
df = rbind(df.1,df.2)

# y1 = log(y1)
# y2 = log(y2)
# y1=y1[is.finite(y1)]
# y2=y2[is.finite(y2)]

p.1 <- ggplot(df.1, aes(x, y)) + 
  geom_point(colour = ifelse(is.finite(df.1$ly),yes = alpha("black",.5),no = alpha("red",.5))) + 
  lims(x = range(c(df$x)), y = range(c(y1,y2))) +
  theme_classic() + 
  theme(axis.ticks = element_line(colour = "black"),
        axis.title = element_text(size = text.size),
        axis.text = element_text(size = text.size, colour = "black"),
        axis.text.x = element_text(size = text.size),
        axis.text.y = element_text(size = text.size))

p.2 <- ggplot(df.2, aes(x, y)) + 
  geom_point(colour = ifelse(is.finite(df.2$ly),yes = alpha("black",.5),no = alpha("red",.5))) + 
  lims(x = range(c(df$x)), y = range(c(y1,y2))) +
  theme_classic() + 
  theme(axis.ticks = element_line(colour = "black"),
        axis.title = element_text(size = text.size),
        axis.text = element_text(size = text.size, colour = "black"),
        axis.text.x = element_text(size = text.size),
        axis.text.y = element_text(size = text.size))
p.1
p.2
# grid.arrange(p.1, p.2, nrow = 1)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
# mean(y1)

p.h1 <- p.1 + geom_hline(yintercept=gm_mean(y1), linetype="dashed", color = col.line, size=size.line)
p.h1
p.h.v1 <- p.h1 + geom_vline(xintercept=mean(x1), linetype="dashed", color = col.line, size=size.line)
p.h.v1
p.h2 <- p.2 + geom_hline(yintercept=gm_mean(y2), linetype="dashed", color = col.line, size=size.line)
p.h2
p.h.v2 <- p.h2 + geom_vline(xintercept=mean(x2), linetype="dashed", color = col.line, size=size.line)
p.h.v2

p.h.v1.s <- p.h.v1 + geom_hline(yintercept=seq(range(y1,y2)[1],range(y1,y2)[2]+1,by = 5), linetype="dotted", color = col.line, size=.3)
p.h.v2.s <- p.h.v2 + geom_hline(yintercept=seq(range(y1,y2)[1],range(y1,y2)[2]+1,by = 5), linetype="dotted", color = col.line, size=.3)
p.h.v1.s
p.h.v2.s


positions.1 <- data.frame(
  x = c(mean(x1)-sd(x1), mean(x1)+sd(x1), mean(x1)+sd(x1), mean(x1)-sd(x1)),
  y = c(mean(y1)-sd(y1), mean(y1)-sd(y1), mean(y1)+sd(y1), mean(y1)+sd(y1))
)
positions.2 <- data.frame(
  x = c(mean(x2)-sd(x2), mean(x2)+sd(x2), mean(x2)+sd(x2), mean(x2)-sd(x2)),
  y = c(mean(y2)-sd(y2), mean(y2)-sd(y2), mean(y2)+sd(y2), mean(y2)+sd(y2))
)

positions.x1 <- data.frame(
  x = c(mean(x1)-sd(x1), mean(x1)+sd(x1), mean(x1)+sd(x1), mean(x1)-sd(x1)),
  y = c(-Inf, -Inf, Inf, Inf)
)
positions.x2 <- data.frame(
  x = c(mean(x2)-sd(x2), mean(x2)+sd(x2), mean(x2)+sd(x2), mean(x2)-sd(x2)),
  y = c(-Inf, -Inf, Inf, Inf)
)
positions.y1 <- data.frame(
  x = c(Inf, -Inf, -Inf, Inf),
  y = c(mean(y1)-sd(y1), mean(y1)-sd(y1), mean(y1)+sd(y1), mean(y1)+sd(y1))
)
positions.y2 <- data.frame(
  x = c(Inf, -Inf, -Inf, Inf),
  y = c(mean(y2)-sd(y2), mean(y2)-sd(y2), mean(y2)+sd(y2), mean(y2)+sd(y2))
)
# p.h.v + geom_polygon(data = positions, aes(x = x, y = y),fill = alpha("gray70",.6))
p.h.v.sd.x1 = p.h.v1.s + geom_polygon(data = positions.x1, aes(x = x, y = y),fill = alpha("gray70",.4))
p.h.v.sd.x1
p.h.v.sd.x2 = p.h.v2.s + geom_polygon(data = positions.x2, aes(x = x, y = y),fill = alpha("gray70",.4))
p.h.v.sd.x2

p.h.v.sd.x.sd.y1 = p.h.v.sd.x1 + geom_polygon(data = positions.y1, aes(x = x, y = y),fill = alpha("gray70",.4))
p.h.v.sd.x.sd.y1
p.h.v.sd.x.sd.y2 = p.h.v.sd.x2 + geom_polygon(data = positions.y2, aes(x = x, y = y),fill = alpha("gray70",.4))
p.h.v.sd.x.sd.y2

p.h.v.sd.x.sd.y.l1 = p.h.v.sd.x.sd.y1 +  geom_smooth(method = "glm",method.args = list(family = "poisson"), se = TRUE, col =alpha("red",.5))
p.h.v.sd.x.sd.y.l1
p.h.v.sd.x.sd.y.l2 = p.h.v.sd.x.sd.y2 +  geom_smooth(method = "glm",method.args = list(family = "poisson"), se = TRUE, col =alpha("red",.5))
p.h.v.sd.x.sd.y.l2

p.marg1 = ggExtra::ggMarginal(p.h.v.sd.x.sd.y.l1, type = "histogram", fill = "gray80", col = "gray70")
p.marg1
p.marg2 = ggExtra::ggMarginal(p.h.v.sd.x.sd.y.l2, type = "histogram", fill = "gray80", col = "gray70")
p.marg2

# grid.arrange(p.marg1,p.marg2, nrow =1)



## ----poissonX_Y1, echo=FALSE, fig.width=5,fig.height=4--------------------------------------------------------------------------------------------------------
p.1


## ----poissonX_Y2, echo=FALSE, fig.width=5,fig.height=4--------------------------------------------------------------------------------------------------------
p.2


## ----poissonX_Y3, echo=FALSE, fig.width=5,fig.height=4--------------------------------------------------------------------------------------------------------
p.marg1


## ----poissonX_Y4, echo=FALSE, fig.width=5,fig.height=4--------------------------------------------------------------------------------------------------------
p.marg2


## ----poissonX_logY, echo=FALSE, fig.show='hide', fig.width=8,fig.height=4-------------------------------------------------------------------------------------
library(scales)
set.seed(1217)
# n = 250
# b0 = 0.6
# b1 = 1.00
# b2 = 0.5
# x1 = rnorm(n, mean = 5, sd = 1) # 
# x1 = runif(n, min = 1, max = 3)
# lambda1 = exp(b0 + b1*x1) # This defines the mean at each point 
# y1 = rpois(n, lambda = lambda1) 
# x2 = rnorm(n, mean = 3, sd = 1)
# lambda2 = exp(b0 + b2*x2)
# y2 = rpois(n, lambda = lambda2) 

y1 = log(y1)
y2 = log(y2)
y1=y1[is.finite(y1)]
y2=y2[is.finite(y2)]

p.1.log <- ggplot(df.1, aes(x, log(y))) + 
  geom_point(colour = ifelse(is.finite(df.1$ly),yes = alpha("black",.5),no = alpha("red",.5))) + 
  lims(x = range(c(df$x)), y = range(c(y1,y2))) +
  theme_classic() + 
  theme(axis.ticks = element_line(colour = "black"),
        axis.title = element_text(size = text.size),
        axis.text = element_text(size = text.size, colour = "black"),
        axis.text.x = element_text(size = text.size),
        axis.text.y = element_text(size = text.size)) #+
  # scale_y_continuous(trans = log_trans(),breaks = c(1:5,seq(10,100,by =10)),limits = c(x = range(c(df$x)), y = range(c(y1,y2))))#,
# breaks = trans_breaks("log", function(x) exp(x)),
# labels = trans_format("log", math_format(e^.x)))


p.2.log <- ggplot(df.2, aes(x, log(y))) + 
  geom_point(colour = ifelse(is.finite(df.2$ly),yes = alpha("black",.5),no = alpha("red",.5))) + 
  lims(x = range(c(df$x)), y = range(c(y1,y2))) +
  theme_classic() + 
  theme(axis.ticks = element_line(colour = "black"),
        axis.title = element_text(size = text.size),
        axis.text = element_text(size = text.size, colour = "black"),
        axis.text.x = element_text(size = text.size),
        axis.text.y = element_text(size = text.size))#+
  # scale_y_continuous(trans = log_trans(),breaks = c(1:5,seq(10,100,by =10)),limits = c(x = range(c(df$x)), y = range(c(y1,y2))))

# grid.arrange(p.1, p.2, nrow = 1)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
# mean(y1)

p.h1 <- p.1.log + geom_hline(yintercept=mean(y1), linetype="dashed", color = col.line, size=size.line)
p.h1
p.h.v1 <- p.h1 + geom_vline(xintercept=mean(x1), linetype="dashed", color = col.line, size=size.line)
p.h.v1
p.h2 <- p.2.log + geom_hline(yintercept=mean(y2), linetype="dashed", color = col.line, size=size.line)
p.h2
p.h.v2 <- p.h2 + geom_vline(xintercept=mean(x2), linetype="dashed", color = col.line, size=size.line)
p.h.v2

p.h.v1.s <- p.h.v1 + geom_hline(yintercept=seq(range(y1,y2)[1],range(y1,y2)[2]+1,by = 5), linetype="dotted", color = col.line, size=.3)
p.h.v2.s <- p.h.v2 + geom_hline(yintercept=seq(range(y1,y2)[1],range(y1,y2)[2]+1,by = 5), linetype="dotted", color = col.line, size=.3)
p.h.v1.s
p.h.v2.s


positions.1 <- data.frame(
  x = c(mean(x1)-sd(x1), mean(x1)+sd(x1), mean(x1)+sd(x1), mean(x1)-sd(x1)),
  y = c(mean(y1)-sd(y1), mean(y1)-sd(y1), mean(y1)+sd(y1), mean(y1)+sd(y1))
)
positions.2 <- data.frame(
  x = c(mean(x2)-sd(x2), mean(x2)+sd(x2), mean(x2)+sd(x2), mean(x2)-sd(x2)),
  y = c(mean(y2)-sd(y2), mean(y2)-sd(y2), mean(y2)+sd(y2), mean(y2)+sd(y2))
)

positions.x1 <- data.frame(
  x = c(mean(x1)-sd(x1), mean(x1)+sd(x1), mean(x1)+sd(x1), mean(x1)-sd(x1)),
  y = c(-Inf, -Inf, Inf, Inf)
)
positions.x2 <- data.frame(
  x = c(mean(x2)-sd(x2), mean(x2)+sd(x2), mean(x2)+sd(x2), mean(x2)-sd(x2)),
  y = c(-Inf, -Inf, Inf, Inf)
)
positions.y1 <- data.frame(
  x = c(Inf, -Inf, -Inf, Inf),
  y = c(mean(y1)-sd(y1), mean(y1)-sd(y1), mean(y1)+sd(y1), mean(y1)+sd(y1))
)
positions.y2 <- data.frame(
  x = c(Inf, -Inf, -Inf, Inf),
  y = c(mean(y2)-sd(y2), mean(y2)-sd(y2), mean(y2)+sd(y2), mean(y2)+sd(y2))
)
# p.h.v + geom_polygon(data = positions, aes(x = x, y = y),fill = alpha("gray70",.6))
p.h.v.sd.x1 = p.h.v1.s + geom_polygon(data = positions.x1, aes(x = x, y = y),fill = alpha("gray70",.4))
p.h.v.sd.x1
p.h.v.sd.x2 = p.h.v2.s + geom_polygon(data = positions.x2, aes(x = x, y = y),fill = alpha("gray70",.4))
p.h.v.sd.x2

p.h.v.sd.x.sd.y1 = p.h.v.sd.x1 + geom_polygon(data = positions.y1, aes(x = x, y = y),fill = alpha("gray70",.4))
p.h.v.sd.x.sd.y1
p.h.v.sd.x.sd.y2 = p.h.v.sd.x2 + geom_polygon(data = positions.y2, aes(x = x, y = y),fill = alpha("gray70",.4))
p.h.v.sd.x.sd.y2

p.h.v.sd.x.sd.y.l1 = p.h.v.sd.x.sd.y1 +  geom_smooth(method = "lm",method.args = list(family = "poisson"), se = TRUE, col =alpha("red",.5))
p.h.v.sd.x.sd.y.l1
p.h.v.sd.x.sd.y.l2 = p.h.v.sd.x.sd.y2 +  geom_smooth(method = "lm",method.args = list(family = "poisson"), se = TRUE, col =alpha("red",.5))
p.h.v.sd.x.sd.y.l2

p.marg1 = ggExtra::ggMarginal(p.h.v.sd.x.sd.y.l1, type = "histogram", fill = "gray80", col = "gray70")
p.marg1
p.marg2 = ggExtra::ggMarginal(p.h.v.sd.x.sd.y.l2, type = "histogram", fill = "gray80", col = "gray70")
p.marg2


## ----poissonX_logY3, echo=FALSE, fig.width=5,fig.height=4-----------------------------------------------------------------------------------------------------
p.marg1


## ----poissonX_logY4, echo=FALSE, fig.width=5,fig.height=4-----------------------------------------------------------------------------------------------------
p.marg2


## ----r_tips1, eval=FALSE--------------------------------------------------------------------------------------------------------------------------------------
## # Description  ------------------------------------------------------------
## # This is the description section as previously presented
## 
## # Libraries ---------------------------------------------------------------
## # Here load the libraries used in the script
## library(ggplot2)
## 
## # Functions ---------------------------------------------------------------
## # Add and describe the functions used in the script
## ## Add more information
## function.name = function(var){tmp=2+var; return(tmp)}
## 
## # Plots -------------------------------------------------------------------
## # Plotting the data simulated
## ## Add more information
## plot(function.name(1:10)+rnorm(10))


## ----set.seed_function, echo=FALSE----------------------------------------------------------------------------------------------------------------------------
set.seed(123)


## ----runif_example--------------------------------------------------------------------------------------------------------------------------------------------
runif(n = 1, min = 1, max = 10) # Gives a random number between 1 and 10
runif(n = 1, min = 1, max = 10) # RNG wasn't reset, different answer (see above)
runif(n = 1, min = 1, max = 10) # Different again... 

set.seed(42); runif(n = 1, min = 1, max = 10) # This sets the RNG 
set.seed(42); runif(n = 1, min = 1, max = 10) # The exact same number 



## ----set.seed_hidden, echo=FALSE------------------------------------------------------------------------------------------------------------------------------
set.seed(123)


## ----sample_numerical_example---------------------------------------------------------------------------------------------------------------------------------
set.seed(12) # Set the RNG 
v.1.10 = 1:10 # Make a vector from 1 to 10 
# Randomly pick 1 (size) value from the vector (x), without replacement 
sample(x = v.1.10, size = 1, replace = FALSE) 


## ----sample_characters_example--------------------------------------------------------------------------------------------------------------------------------
set.seed(3) # Set the RNG 
# Randomly pick 5 (size) letters from the vector (x), without replacement 
sample(x = LETTERS, size = 5, replace = FALSE) 
sample(x = as.factor(month.abb), size = 5, replace = FALSE) 


## ----permutations_load_viridis, echo=FALSE--------------------------------------------------------------------------------------------------------------------
library(viridis)


## ----permutations_df, fig.width=4,fig.height=3----------------------------------------------------------------------------------------------------------------
set.seed(123)
n = 40; col = viridis::viridis(n = n)
x = 1:n ; y = 2+.5*x + rnorm(n, sd=7)
df.xy = data.frame(x,y, col )


## ----permutations_XY, fig.width=4,fig.height=3----------------------------------------------------------------------------------------------------------------
set.seed(321)
df.xy$x.s = sample(df.xy$x) #<<
df.xy$y.s = sample(df.xy$y) #<<
# We break up the link of X and Y 


## ----permutations_plot, echo=-c(1:5), fig.width=12,fig.height=4-----------------------------------------------------------------------------------------------
par(mfrow=c(1,3), mar=c(3,4,2,1), cex = 1.2)
plot.lm.c <- function(data, formula,main = "") {
  plot(as.formula(formula), col=col, xlab ="", data=data, pch=19, main = main)
  abline(lm(as.formula(formula),  data=data))
  title(xlab = as.character(as.formula(formula)[3]), line=2, cex.lab=1.2)}
plot.lm.c(df.xy,"y~x",    main = "Original x-y")
plot.lm.c(df.xy,"y~x.s",  main = "Permutated x")
plot.lm.c(df.xy,"y.s~x.s",main = "Permutated x and y")


## ----perm.boot, eval=FALSE------------------------------------------------------------------------------------------------------------------------------------
## permutate.df = replicate(n = 200, # nperm
##           expr = data.frame(y = sample(df.lm$y,size = nrow(df.lm), replace = FALSE), #<<
##                             x = df.lm$x), simplify = FALSE)
## bootstrap.fun <- function(data) {
##   tmp.sample = sample(1:nrow(data),size = c(nrow(data)-20), replace = FALSE) #<<
##   data.frame(y = data[tmp.sample,"y"], x = data[tmp.sample,"x"]) } # end boot fun
## bootstrap.df = replicate(n = 200, expr = bootstrap.fun(df.lm), simplify = FALSE)


## ----permutations_bootstrap_df_lm, echo=FALSE, results='hide', message=FALSE, warning=FALSE, fig.width=7, fig.height=3.5--------------------------------------
par(mar=c(4,4,0.1,0.1))
set.seed(12345678)
n = 100; beta0 = 2.5; beta1 = 0.8
x.lm = rnorm(n = n, mean = 10, sd = 1)
err = rnorm(n = n, mean = 0, sd = 1)
# Linear combination 
y.lm = beta0 + beta1*x.lm + err
# Make a dataframe of the data 
df.lm = data.frame(x = x.lm, y = y.lm)
par(mar = c(4,4,.5,.5))
# Colour 
b.5 = scales::alpha("black",alpha = .5)

# PLot the data 
plot(y~x, data = df.lm, pch = 19, col = b.5, xlab = "")
title(xlab = "x", ylab="", line=2.2, cex.lab=1.2)
# Model the data 
lm.out = lm(y~x, data = df.lm)
# Add a line to the plot 
abline(lm.out)
pred = predict(lm.out, interval = "confidence") # 95%  confidence interval
polygon(x = c(df.lm$x[order(df.lm$x)], rev(df.lm$x[order(df.lm$x)])),  
        y = c(pred[, "lwr"][order(df.lm$x)],rev(pred[, "upr"][order(df.lm$x)])), 
        col = scales::alpha("blue",.5))



sim.lm = simulate(lm.out, nsim = 1e5, seed = 12) %>% data.frame
lower_ci_sim <- apply(sim.lm, 1, function(x) quantile(x, probs = 0.025) )
upper_ci_sim <- apply(sim.lm, 1, function(x) quantile(x, probs = 0.975) )
sims_summary <- data.frame(
  lower = lower_ci_sim,
  upper = upper_ci_sim
)


# Add permutations
set.seed(98765432)
permutate.df = replicate(n = 200, # nperm 
          expr = data.frame(y = sample(df.lm$y,size = nrow(df.lm), replace = FALSE), x = df.lm$x),
          simplify = FALSE)
lm.out.perm = mapply(lm, permutate.df)
apply(lm.out.perm, 2,function(x) abline(x, col = scales::alpha("orange",.5)))

# Add bootstrap
bootstrap.fun <- function(data) {
  tmp.sample = sample(1:nrow(data),size = c(nrow(data)-20),replace = FALSE)
  data.frame(y = data[tmp.sample,"y"], 
             x = data[tmp.sample,"x"])
}
set.seed(54)
bootstrap.df = replicate(n = 200, # nperm 
                        expr = bootstrap.fun(df.lm),
                        simplify = FALSE)

lm.out.boot = mapply(lm, bootstrap.df)
apply(lm.out.boot, 2, function(x) abline(x, col = scales::alpha("red",.9)))

abline(lm.out)
points(df.lm$x,df.lm$y, col = b.5, pch = 19)
points(mean(df.lm$x),mean(df.lm$y), col = "green", pch = 19)
abline(v =mean(df.lm$x), h = mean(df.lm$y), lty = 3)



## ----rdm_dates------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(42)
# Get a sequence of dates
datae.seq = seq(from = as.Date('2010/01/01'), 
                to = as.Date('2022/01/01'), 
                by = "day")

# Look at beginning and end of sequence 
head(datae.seq, 4); tail(datae.seq, 4)

# Get only 5 elements of the generated sequence 
sample(datae.seq, 5)


## ----sample_die, echo=-1, fig.width=8, fig.height=3.5---------------------------------------------------------------------------------------------------------
par(mar=c(4,4,1,0.1))
set.seed(50)
p_dice = c(1,1,1,1,1,5) # Here we have twice the change of landing on 6
                        # Same as writing p_dice/sum(p_dice) or the prob.
nb.tosses = 100
die_results <- sample(x = 1:6, # or seq(from = 1, to=6, by=1)
                      size = nb.tosses,
                      replace = T, prob = p_dice) 
barplot(table(die_results), ylab = "Fq", xlab ="Face", main = "Loaded dice") # table(die_results)/nb.tosses


## ----rep_function_example-------------------------------------------------------------------------------------------------------------------------------------
(let4first = LETTERS[1:4])
rep(let4first, times = 2) # Repeat the WHOLE sequence twice 
rep(let4first, each = 2) # Repeat each element twice 

# Set the number of repeat for each element in the vector 
rep(let4first, times = c(1,2,3,6))

# Complete replication: replicate each element twice and do this three times 
rep(let4first, each = 2, times = 3)
rep(let4first, length.out = 6) # Repeat the vector until you hit the length.out



## ----gl_func--------------------------------------------------------------------------------------------------------------------------------------------------
nb.of.levels = 2
nb.of.replicates = 8
labels.for.the.factors = c("Control", "Treat")

## First control, then treatment:
gl(n = nb.of.levels, k = nb.of.replicates, labels = labels.for.the.factors)

## 20 alternating As and Bs
gl(n = 2, k = 1, length = 20, labels = LETTERS[1:2])

## alternating pairs of 1s and 2s
gl(n = 2, k = 2, length = 19) # see last element missing


## ----data_replicate-------------------------------------------------------------------------------------------------------------------------------------------
set.seed(1234)
data.replicated = replicate(n = 2,
                            expr = data.frame(gr = rep(LETTERS[1:3], each = 2),
                                              y = rnorm(6)), 
                            simplify = FALSE)


## ----data_replicate1------------------------------------------------------------------------------------------------------------------------------------------
data.replicated[[1]]


## ----data_replicate2------------------------------------------------------------------------------------------------------------------------------------------
data.replicated[[2]]


## ----expand.example-------------------------------------------------------------------------------------------------------------------------------------------
set.seed(1234)
exp.design = expand.grid(height = seq(60, 70, 10), 
                         weight = seq(100, 200, 100),
                         sex = c("Male","Female"), 
                         stringsAsFactors = TRUE) # characters will be factors
exp.design


## ----deck_of_caRds, echo=-1-----------------------------------------------------------------------------------------------------------------------------------
# Nicer deck 
cards = c( "A",2:10, "J", "Q", "K"); length(cards) # Get cards type
suits = c("♦", "♣", "♥", "♠" ) # Get the suits c("Diamonds", "Clubs", "Hearts", "Spades")
cute.deck <- expand.grid(cards = cards, suits = suits) #<<
cute.deck$nice.cards = apply(cute.deck, 1, function(x) paste0(x, collapse = "")) # Combine cards and suits
cute.deck$col = ifelse(cute.deck$suits %in% c("♦","♥"),"red","black") # add colour 

# Select cards at random 
set.seed(1234)
n.c = 5
row.sel = sample(1:nrow(cute.deck),size = n.c, replace = FALSE) #<<
cute.deck[row.sel,"nice.cards"]


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
# probability of extracting 2 aces from a deck of 52 cards which has 4 aces and that we draw 3 cards
nb.aces.in.deck = 4
nb.of.cards.in.deck = 52
draw.cards = 3
nb.cards.wish.to.check = 2
# This is modeled with a hypergeometric distribution (discrete), which calculate probabilities when sampling without replacement
dhyper(nb.cards.wish.to.check, 
       m = nb.aces.in.deck,
       n = nb.of.cards.in.deck-nb.aces.in.deck,
       k = draw.cards)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Check prop. of each suits in a deck 
table(cute.deck$suits)/nrow(cute.deck) 

# cards <- c("Ace", "Deuce", "Three", "Four","Five", "Six", "Seven", "Eight", "Nine", "Ten", "Jack", "Queen", "King")
# Define suits, cards, values
suits <- c("Diamonds", "Clubs", "Hearts", "Spades")
cards <- c("Ace", 2:10, "Jack", "Queen", "King")
# Build deck, replicated proper number of times
deck <- expand.grid(cards=cards, suits=suits)
deck$value <- c(1, 2:9, rep(10, 4))
deck$col = ifelse(deck$suits %in% c("Diamonds","Hearts"),"red","black")


## ----outer.example.nt, eval=FALSE-----------------------------------------------------------------------------------------------------------------------------
## x <- 1:10; names(x) <- x
## # Multiplication & Power Tables
## x %o% x # same as outer(x, x, "*")
## x.mul = outer(x, x, "*")
## x.sub = outer(x, x, "-")
## x.add = outer(x, x, "+")
## x.div = outer(x, x, "/")
## 


## ----outer.example, echo=FALSE, fig.width=8,fig.height=5------------------------------------------------------------------------------------------------------
x <- 1:10; names(x) <- x; y <- 1:10; names(y) <- y
# Multiplication & Power Tables
# x %o% x # same as outer(x, x, "*")
x.mul = outer(x, y, "*")
x.sub = outer(x, y, "-")
x.add = outer(x, y, "+")
x.div = outer(x, y, "/")
add.lab = function(x, out.mat, cex = 1) {
  for (i in 1:length(x)) {
    for (j in 1:length(x)) {
      text(i,j,labels = round(out.mat[i,j],2),cex = cex)
    }
  }
}
cex=.8
par(mfrow = c(2,2), mar = c(4,4,1,1))
image(x = x, y = y, z = x.mul, axes = FALSE, main ='Multiplication', xlab = "")
title(xlab = "x", line=2, cex.lab=1.2)
axis(1, at = seq(1, max(x), by = 1))
axis(2, at = seq(1, max(y), by = 1))
add.lab(x,x.mul, cex=cex)

image(x = x, y = y, z =x.sub, axes = FALSE, main ='Subtraction', xlab = "")
title(xlab = "x", line=2, cex.lab=1.2)
axis(1, at = seq(1, max(x), by = 1))
axis(2, at = seq(1, max(y), by = 1))
add.lab(x,x.sub, cex=cex)

image(x = x, y = y, z =x.add, axes = FALSE, main ='Addition', xlab = "")
title(xlab = "x", line=2, cex.lab=1.2)
axis(1, at = seq(1, max(x), by = 1))
axis(2, at = seq(1, max(y), by = 1))
add.lab(x,x.add, cex=cex)

image(x = x, y = y, z =x.div, axes = FALSE, main ='Division', xlab = "")
title(xlab = "x", line=2, cex.lab=1.2)
axis(1, at = seq(1, max(x), by = 1))
axis(2, at = seq(1, max(y), by = 1))
add.lab(x,x.div, cex=cex)

y <- 2:8; names(y) <- paste(y,":", sep = "")
# outer(y, x, "^")

# outer(month.abb, 1999:2003, FUN = "paste")

## three way multiplication table:
# x %o% x %o% y[1:3]



## ----wallpapeR, eval=FALSE------------------------------------------------------------------------------------------------------------------------------------
## #<<
## par(bg = NA, mar = c(0,0,0,0)) # will take the whole screen
## x = 1:70; y = 1:70 # make variables
## x.val = (x^2);   y.val = y
## xy.out = outer(x.val, y.val, "-") # manipulate the variables
## nb.col.pal = length(grDevices::hcl.pals()) # find all palettes in hcl.pals
## #<<
## 
## xy.out = xy.out + rnorm(length(xy.out), mean = 0, sd = 200) # Prettify with RDM
## image(t(xy.out),xaxt = "n", yaxt = "n", bty = "n",
##       col = hcl.colors(1000, palette = hcl.pals()[rdm.nb])) # 61,93 looks good


## ----waller, echo=FALSE, eval=FALSE---------------------------------------------------------------------------------------------------------------------------
## # Make computer wallpapers
## # png("~/Desktop/wallpapers_sunwaves.png", res = 300, units = "px",
## #     width = 1920, # add your screen resolution here
## #     height = 1080)
## set.seed(9)
## par(bg = NA, mar = c(0,0,0,0)) # will take the whole screen
## x = 1:70; y = 1:70 # make variables
## x.val = (x^2);   y.val = y
## xy.out = outer(x.val, y.val, "-") # manipulate the variables
## nb.col.pal = length(grDevices::hcl.pals()) # find all palettes in hcl.pals
## rdm.nb = sample(c(1:nb.col.pal), size = 1); print(rdm.nb) # Get random number
## xy.out = xy.out + rnorm(length(xy.out), mean = 0, sd = 200) # Prettify with RDM
## # f.na.iamge = which(is.na(xy.out))
## # xy.out[f.na.iamge] <- rnorm(length(f.na.iamge),mean = mean(tan(x)))
## image(t(xy.out),
##       col = hcl.colors(1000, palette = hcl.pals()[rdm.nb]), # 61,93 looks good
##       xaxt = "n", yaxt = "n", bty = "n")
## # dev.off()


## ----wallpapeR_solution, fig.width=6, fig.height=2------------------------------------------------------------------------------------------------------------
set.seed(9) #<<
par(bg = NA, mar = c(0,0,0,0)) # will take the whole screen
x = 1:70; y = 1:70 # make variables 
x.val = (x^2);   y.val = y
xy.out = outer(x.val, y.val, "-") # manipulate the variables
nb.col.pal = length(grDevices::hcl.pals()) # find all palettes in hcl.pals
rdm.nb = sample(c(1:nb.col.pal), size = 1); print(rdm.nb) # Get random number #<<
xy.out = xy.out + rnorm(length(xy.out), mean = 0, sd = 200) # Prettify with RDM
image(t(xy.out),xaxt = "n", yaxt = "n", bty = "n",
      col = hcl.colors(1000, palette = hcl.pals()[rdm.nb])) # 61,93 looks good


## ----show_sample, eval=FALSE----------------------------------------------------------------------------------------------------------------------------------
## sample()


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(12345) # Sets the random number generator to a fix value
vec.1.10 = 1:10 # Make the vector to choose from 
mean(vec.1.10)
size = 10
n = 10
replicate.draws = replicate(n = n,simplify = TRUE,
          expr = sample(x = vec.1.10, size = size, replace = TRUE))

get.average = apply(replicate.draws, 2, mean)
hist(get.average)

# Check out what simplify does and deal with it 
replicate.draws.list = replicate(n = n,simplify = FALSE,
                            expr = sample(x = vec.1.10, size = size, replace = TRUE))
# Make into a matrix
do.call(cbind, replicate.draws.list)
do.call(rbind, replicate.draws.list)


## ----sample_replacement---------------------------------------------------------------------------------------------------------------------------------------
set.seed(12345) # Sets the random number generator to a fix value
vec.1.10 = 1:10 # Make the vector to choose from 
sample(x = vec.1.10, size = 4, replace = FALSE) # Sample 4 nb without replacement
sample(x = vec.1.10, size = 4, replace = TRUE) # Sample 4 nb with replacement


## ----sample_example-------------------------------------------------------------------------------------------------------------------------------------------
set.seed(123); table(sample(x = 0:1, size = 1000, replace = T,prob = c(.5,.5)))


## ----linear_no_error, echo=-1, fig.width=6, fig.height=5------------------------------------------------------------------------------------------------------
par(mar=c(4,4,.5,.5), cex = 1.5)
x = 1:10
y = 2 + 3 * x
gr = rep(letters[1:2],each = 5)
linear.df = data.frame(x,y,gr)
plot(y~x, col = as.factor(gr), data = linear.df, pch = 19, cex = 1.5)


## ----coin_tosses, echo=FALSE, fig.width=10,fig.height=4-------------------------------------------------------------------------------------------------------
# This is an example of "long run frequency" probability
par(mfrow = c(1,2), mar = c(4,4,2.5,.5), cex = 1.2)
toss.coin <- function(n.tosses = 1000) {
  ctos = rbinom(n.tosses, 1, .5)
  nb = 1:length(ctos)
  r.sum = cumsum(ctos)
  d.coin.toss = data.frame(nb, ctos, r.sum, run.prop = r.sum/nb)
  d.coin.toss
}
set.seed(123456)
rep.toss = replicate(n = 3,toss.coin(n.tosses = 10), simplify = F)
result = lapply(rep.toss, "[", "run.prop")
tosses= do.call(cbind,result)
matplot(tosses,type = "l", ylim = c(0,1), main = "3 coins, 10 tosses",
        ylab = "Proportion of tosses", xlab = "Number of tosses")
abline(h = 0.5,lty = 3, lwd= 3)

set.seed(123456)
rep.toss = replicate(n = 100,toss.coin(n.tosses = 5e3), simplify = F)
result = lapply(rep.toss, "[", "run.prop")
tosses= do.call(cbind,result)
matplot(tosses,type = "l", ylim = c(0,1),  main = "100 coins, 5000 tosses",
        ylab = "Proportion of tosses", xlab = "Number of tosses")
abline(h = 0.5,lty = 3, lwd= 3)

# plot(d.coin.toss$run.prop~d.coin.toss$nb, type = "l", ylim = c(0,1))


## ----odds_prob_scale_zero_one, echo=FALSE, fig.width=4.7,fig.height=2-----------------------------------------------------------------------------------------
par(mar = c(0,0,0,0), 
    # mfrow=c(1,2),
    cex = 1.2)
# Probability
plot(0, type ="n", xlim = c(0,1), axes = F,ylab = "", xlab = "")
segments(0,0,1,0)
yval = .1
lwd = 2
seqp=seq(0,1,by = .1)
for (i in seqp) {
  segments(i,-yval,i,yval)
}
ytext = .65
text(seqp,y = -.25,labels = seqp)
text(.5,y = ytext,labels = c("Probability"))
arrows(0,.3,1,.3,length = .1,code = 3, lwd = lwd, col = "gray10")



## ----odds_scale_zero_inf, echo=FALSE, fig.width=4.7,fig.height=2----------------------------------------------------------------------------------------------
par(mar = c(0,0,0,0), 
    # mfrow=c(1,2),
    cex = 1.2)
# Odds
plot(0, type ="n", xlim = c(0,11), axes = F,ylab = "", xlab = "")
segments(0,0,10,0)
yval = .1
lwd = 2
for (i in 0:10) {
  segments(i,-yval,i,yval)
}
text(0:10,y = -.25,labels = c(0:10))
text(.5,y = ytext,labels = c("Odds \nlosing"))
text((10+1)/2,y = ytext,labels = c("Odds \nwinning"))
arrows(0,.3,1,.3,length = .1,code = 3, lwd = lwd, col = "red")
arrows(1,.3,11,.3,length = .1,code = 3, lwd = lwd, col = "blue")


## ----log_odds_scale_inf_inf, echo=FALSE, fig.width=4.7,fig.height=2-------------------------------------------------------------------------------------------
par(mar = c(0,0,0,0), cex = 1.2)
plot(0, type ="n", xlim = c(-6,6), axes = F,ylab = "", xlab = "")
segments(-5,0,5,0)
yval = .1
lwd = 2
for (i in -5:5) {
  segments(i,-yval,i,yval)
}
ytext = .65
text(-5:5,y = -.25,labels = c(-5:5))
text(-5/2,y = ytext,labels = c("Log odds \nlosing"))
text((5)/2,y = ytext,labels = c("Log odds \nwinning"))
arrows(-6,.3,0,.3,length = .1,code = 1, lwd = lwd, col = "red")
arrows(0,.3,6,.3,length = .1,code = 2, lwd = lwd, col = "blue")
points(0,.3, pch =19)



## ----prob_odd_log_trans---------------------------------------------------------------------------------------------------------------------------------------
raw.num = c(1,2,8) # numbers NOT on the log scale (simulate a DOUBLING process like qPCR)
mean(raw.num) # Arithmetic average (which would be WRONG)
mean.log = mean(log(raw.num)/log(2)) # Taking the mean of the log (this is the geometric mean)
mean.log # Mean, on the log scale, for the raw numbers transformed: log(raw.num)
2^mean.log
# gm = exp(mean(log(raw.num))) # calculating the geometric mean 
gm = 2^(mean(log(raw.num)/log(2))) # Here we calculate the geometric mean with a base 2
# Now compare the geometric mean of the raw numbers with the mean of the log, transformed back into the raw number scale 
gm
2^mean.log


## ----odds_example, echo=-1------------------------------------------------------------------------------------------------------------------------------------
library(MASS)
tot.nb.ev = 100; success = 20
failure = tot.nb.ev-success
p = success/tot.nb.ev; q = failure/tot.nb.ev
(odds.favor = (success)/(failure)) # Odds in favor of event ( 1 in 4), same as p/q
(odds.agnst = (failure)/(success)) # Odds against the event (here 4 to 1), same as q/p


## ----coin_sequence_image, echo=FALSE, fig.width=4.7,fig.height=2----------------------------------------------------------------------------------------------
par(mar = c(0,0,0,0))
# Get circles
nb.flips = 5
radius =1 
# initialize a plot
plot(x = c(0, nb.flips*(radius+.5)+5), y = c(-1, 1), type = "n", 
     axes = F,
     asp = 1,ylab = "", xlab = "")
w = 0
pos=0
col.v = c("white",rep("black",nb.flips-1))
for (i in 1:nb.flips) {
  # prepare "circle data"
  radius = 1
  center_x = w + 1
  center_y = pos
  theta = seq(0, 2 * pi, length = 200) # angles for drawing points around the circle
  # draw the circle
  lines(x = radius * cos(theta) + center_x, y = radius * sin(theta) + center_y)
  polygon(x = radius * cos(theta) + center_x, y = radius * sin(theta) + center_y, col = col.v[i])
  text(center_x,center_y+1.5,labels = c("Win",rep("Lose",nb.flips-1))[i], cex = 2)
  w = w + 2*radius+.5
}



## ----probability_example--------------------------------------------------------------------------------------------------------------------------------------
tot.nb.ev = 100; success = 20
failure = tot.nb.ev-success
(probability.favor = (success)/(tot.nb.ev)) # probability of event (here 20%), fractions(probability.favor)
(q = failure/tot.nb.ev)


## ----get_head_tail, echo=FALSE--------------------------------------------------------------------------------------------------------------------------------
ht <- function(d, m=5, n=m){
  # print the head and tail together
  get.ht = list(head = head(d,m), tail = tail(d,n))
  do.call(rbind,get.ht)
}


## ----odds_prob, fig.width=9,fig.height=7.5, highlight.output=c(TRUE, FALSE)-----------------------------------------------------------------------------------
# Get some LOG ODDS numbers 
log_odds = seq(from = -5, to = 5, by = .25)
odds = exp(log_odds) # Transformed into odds 

# Make inverse logit function 
inv.logit = function(x) {exp(x)/(1 + exp(x))} # takes log odds as input. Same as function(x){1/(1 + exp(-x))}

p = inv.logit(log_odds) # This is  p = odds/(1 + odds) 
q = 1-p # Probability of failure (1-p)

# Store log_odds other data to plot 
d = data.frame(log_odds, odds, p, q) 
format(x = ht(d,3),digits = 2, scientific = F)


## ----odds_prob_plot_with_pol, echo=FALSE, eval=TRUE, fig.width=12,fig.height=6--------------------------------------------------------------------------------
# adde prob scale
par(mfrow = c(2,3))
par(mar = c(0,4,0,0), 
    # mfrow=c(1,2),
    cex = 1.2)
# Probability
plot(0, type ="n", xlim = c(0,1), axes = F,ylab = "", xlab = "")
segments(0,0,1,0)
yval = .05
lwd = 2
axis.t = -.17
seqp=seq(0,1,by = .2)
for (i in seqp) {
  segments(i,-yval,i,yval)
}
ytext = .55
text(seqp,y = axis.t,labels = seqp)
text(.5,y = ytext,labels = c("Probability"))
arrows(0,.3,1,.3,length = .1,code = 3, lwd = lwd, col = "gray10")
###
# Odds
plot(0, type ="n", xlim = c(0,11), axes = F,ylab = "", xlab = "")
segments(0,0,10,0)
lwd = 2
for (i in 0:10) {
  segments(i,-yval,i,yval)
}
text(0:10,y = axis.t,labels = c(0:10))
text(.7,y = ytext,labels = c("Odds \nlosing"))
text((10+1)/2,y = ytext,labels = c("Odds \nwinning"))
arrows(0,.3,1,.3,length = .1,code = 3, lwd = lwd, col = "red")
arrows(1,.3,11,.3,length = .1,code = 3, lwd = lwd, col = "blue")
###
plot(0, type ="n", xlim = c(-6,6), axes = F,ylab = "", xlab = "")
segments(-5,0,5,0)
lwd = 2
for (i in -5:5) {
  segments(i,-yval,i,yval)
}
text(-5:5,y = axis.t,labels = c(-5:5))
text(-5/2,y = ytext,labels = c("Log odds \nlosing"))
text((5)/2,y = ytext,labels = c("Log odds \nwinning"))
arrows(-6,.3,0,.3,length = .1,code = 1, lwd = lwd, col = "red")
arrows(0,.3,6,.3,length = .1,code = 2, lwd = lwd, col = "blue")
points(0,.3, pch =19)

###

par(mar =c(4,4,1.5,.25), cex = 1.3)
p.o = d$p~d$odds
p.lo = d$p~d$log_odds
o.lo = d$odds~d$log_odds
plot(o.lo, type="l", ylab="Odds", xlab="", lwd=3, main = "ln(odds) ~ odds"); abline(v=0, lty=3)
polygon(x = c(-10,0,0,-10),y = c(0,0,1,1), col = scales::alpha("red",.5), border =  scales::alpha("red",.5))
polygon(x = c(0,10,10,0),y = c(1,1,150,150), col = scales::alpha("blue",.5), border =  scales::alpha("blue",.5))
title(xlab = "Ln odds", ylab="", line=2.2, cex.lab=1.2)

plot(p.o,  type="l", ylab="Probability", xlab="",    lwd=3, main = "p ~ odds"); abline(h=.5, v=0, lty=3)
polygon(x = c(0,1,1,0),y = c(0,0,.5,.5), col = scales::alpha("red",.5), border =  scales::alpha("red",.5))
polygon(x = c(1,150,150,1),y = c(0.5,0.5,1,1), col = scales::alpha("blue",.5), border =  scales::alpha("blue",.5))
title(xlab = "Odds", ylab="", line=2.2, cex.lab=1.2)

plot(p.lo, type="l", ylab="Probability", xlab="", lwd=3, main = "p ~ ln(odds)"); abline(h=.5, v=0, lty=3)
polygon(x = c(-10,0,0,-10),y = c(0,0,.5,.5), col = scales::alpha("red",.5), border =  scales::alpha("red",.5))
polygon(x = c(0,10,10,0),y = c(0.5,0.5,1,1), col = scales::alpha("blue",.5), border =  scales::alpha("blue",.5))
title(xlab = "Ln odds", ylab="", line=2.2, cex.lab=1.2)


## ----probability_normal_pval, echo=FALSE, eval=TRUE, fig.width=9,fig.height=3---------------------------------------------------------------------------------
par(mfrow = c(1,2), mar=c(3,3,0.1,0.1))
draw.normal(where = "both",mean = 110,  prob = 0.025,ylim = c(0,.6), col = scales::alpha(c('blue',NA),.5), text = TRUE, text.height = .56)
draw.normal(where = "middle2",mean = 110,  prob = 0.4998,ylim = c(0,.6), col = scales::alpha(c('blue',NA),.5), text = F, text.height = .56)


## ----feel_random_var, echo=FALSE------------------------------------------------------------------------------------------------------------------------------
set.seed(123)
n = 1000
rnorm.val = rnorm(n = n, mean = 15, sd = 2)
random.var.normal = round(rnorm.val,1)


## ----get_rdm_nb_1, echo=FALSE---------------------------------------------------------------------------------------------------------------------------------
random.var.normal[1]


## ----get_rdm_nb_5, echo=FALSE---------------------------------------------------------------------------------------------------------------------------------
random.var.normal[1:5]


## ----feel_random_var_hist_5, echo=FALSE, fig.width=6, fig.height=5--------------------------------------------------------------------------------------------
par(mfrow=c(1,1),cex = 1.5)
hist(random.var.normal[1:5], main = "Hisogram of random variable", xlab = "x")


## ----get_rdm_nb_100, echo=FALSE-------------------------------------------------------------------------------------------------------------------------------
random.var.normal[1:100]


## ----feel_random_var_hist_100, echo=FALSE, fig.width=6, fig.height=5------------------------------------------------------------------------------------------
par(mfrow=c(1,1),cex = 1.5)
hist(random.var.normal[1:100], main = "Hisogram of random variable", xlab = "x")


## ----happy_birthday, echo=c(-1), eval=TRUE--------------------------------------------------------------------------------------------------------------------
# from the ??birthday page in R help 
# probability of > 2 people with the same birthday
(a = pbirthday(23, coincident = 3))
# 0.9 probability of >=3 coincident birthdays
(a = qbirthday(coincident = 3, prob = 0.9)) # Gives the number of people
# Chance of >=4  coincident birthdays in 300 people
(a = pbirthday(300, coincident = 4))


## ----pbirthday_plot1------------------------------------------------------------------------------------------------------------------------------------------
## from Diaconis & Mosteller p. 858. 'coincidence' is that husband, wife, daughter all born on the 16th
qbirthday(classes = 30, coincident = 3) # approximately 18

coin = 2
n.seq = 1:qbirthday(p = 0.999,coincident = coin)
get.probs = mapply(pbirthday, n = n.seq, classes = 365, coincident = coin)
plot(get.probs~n.seq, type = "l"); abline(h = 0.5, v = qbirthday(p = 0.5,coincident = coin), lty = 3)


## ----pbirthday_plot2_derivation-------------------------------------------------------------------------------------------------------------------------------
# birthday paradox
# scaling exponents in mathematics
n.people = 23
day.of.year = sample(1:365,n.people, replace = TRUE)
day.of.year
happy.birthday = function(n.people = 20) {
  day.of.year = sample(1:365,n.people, replace = TRUE)
  length(which(table(day.of.year)>2))
}
barplot(table(replicate(10,happy.birthday(n.people = 40))))
pbirthday(40,classes = 365)
prod(c(1-(1:182/365)))

hbirthd <- function(n.p = 5) {
  1-prod(365:c(365-c(n.p-1)))/365^n.p
}
n.p = 2:65
plot(x = n.p,mapply(hbirthd, n.p = n.p), type = "l");abline(h = 0.5, v = qbirthday(), lty = 3)


## ----Statistical_dist, echo=FALSE, fig.width=14,fig.height=8--------------------------------------------------------------------------------------------------
par(mfrow = c(3,4), 
    # bg = NA,
    mar= c(3.5,4,3,3))
# par(mfrow = c(1,1), bg = NA)
col = scales::alpha("black",.5)
col.r = scales::alpha("red",.5)
col.g = scales::alpha("green",.5)
col.b = scales::alpha("blue",.5)
lwd=2
# Binom
#define range of "successes"
success <- 0:20
plot(success, dbinom(x = success, size=1, prob=0.5),type='h', col = col.b, lwd=lwd, ylim = c(0,1), xlab ="",ylab = "Probability", main = "Binomial")
curve(expr = dbinom(x = x, size=20, prob=0.6),col = col, lwd=lwd, ylim = c(0,1), type = 'h',add = T)
curve(expr = dbinom(x = x, size=20, prob=0.4),col = col.r, lwd=lwd, ylim = c(0,1), type = 'h',add = T)
abline(h=seq(0,1, by = .1), lty = 3, lwd = .3)
title(xlab = "Success",line=2.2, cex.lab=1.2) 
legend("topright",legend = c("n = 1, p = .5","n = 20, p = .6","n = 20, p = .4"),lty = c(1,1,1), col =c(col.b,col, col.r), lwd = 2)

# Poisson
plot(success, dpois(success, lambda=1), type='h', col = col.b, lwd=lwd, ylim = c(0,1), xlab = "", ylab = "Probability", main = "Poisson")
points(success, dpois(success, lambda=5), type='h', col=col, lwd=lwd)
points(success, dpois(success, lambda=12.56), type='h', col=col.r, lwd=lwd)
abline(h=seq(0,1, by = .1), lty = 3, lwd = .3)
title(xlab = "Success",line=2.2, cex.lab=1.2) 
legend("topright",legend = c("lambda = 1","lambda = 5","lambda = 12.56"),lty = c(1,1,1), col =c(col.b,col, col.r), lwd = 2)

# Uniform discrete
# Define uniform discrete 
dunifdisc<-function(x, min=0, max=1) ifelse(x>=min & x<=max & round(x)==x, 1/(max-min+1), 0)
curve(dunifdisc(x, 7,10), type = "h",from=6, to=11, col=col, xlab = "",
      lwd=3, ylim = c(0,1), ylab = "Probability", main = "Uniform discrete"); title(xlab = "x",line=2.2, cex.lab=1.2) 
curve(dunifdisc(x, 6,7), type = "h",add=TRUE, col=col.b, xlab = "",
      lwd=2, ylim = c(0,1), ylab = "Probability", main = "Uniform discrete"); title(xlab = "x",line=2.2, cex.lab=1.2) 
abline(h=seq(0,1, by = .1), lty = 3, lwd = .3)
legend("topright",legend = c("min = 6, max = 7", "min = 7, max = 10"),lty = c(1), col =c(col.b, col), lwd = 2)

# Hypergeometric
xs <- 0:20
plot(xs,dhyper(xs, m =20,n = 10,k = 10), type = "h",from=0, to=15, col=col.b, xlab = "",
      lwd=1, ylim = c(0,0.5), ylab = "Probability", main = "Hypergeometric")
curve(expr = dhyper(x, m =30,n = 30,k = 30),col = col, lwd=lwd, ylim = c(0,1), type = 'h',add = T)
title(xlab = "x",line=2.2, cex.lab=1.2) 
abline(h=seq(0,1, by = .1), lty = 3, lwd = .3)
legend("topright",legend = c("m = 20, n = 10, k = 10", "m = 30, n = 30, k = 30"),lty = c(1), col =c(col.b, col), lwd = 2)



# Chi-sq
curve(dchisq(x, df = 1), from = 0, to = 40, col = col.b, lwd=lwd, ylim = c(0,1), xlab = "", ylab = "Density", main = "Chi-square")
curve(dchisq(x, df = 4), from = 0, to = 40, col = col, lwd=lwd, ylim = c(0,1), add = T)
curve(dchisq(x, df = 10), from = 0, to = 40, col = col.r, lwd=lwd, ylim = c(0,1), add = T)
abline(h=seq(0,1, by = .1), lty = 3, lwd = .3)
title(xlab = "x",line=2.2, cex.lab=1.2) 
legend("topright",legend = c("df = 1","df = 4","df = 10"),lty = c(1,1,1), col =c(col.b,col, col.r), lwd = 2)

# Exponential 
curve(dexp(x, rate = .5), from=0, to=10, col=col.b, lwd=lwd, ylim = c(0,1), xlab = "",ylab = "Density", main = "Exponential")
curve(dexp(x, rate = .2), from=0, to=10, col=col, lwd=lwd, ylim = c(0,1),add = T)
curve(dexp(x, rate = .8), from=0, to=10, col=col.r, lwd=lwd, ylim = c(0,1),add = T)
abline(h=seq(0,1, by = .1), lty = 3, lwd = .3)
title(xlab = "x",line=2.2, cex.lab=1.2) 
legend("topright",legend = c("rate = 0.8","rate = 0.5","rate = 0.2"),lty = c(1,1,1), col =c(col.b,col,col.r), lwd = 2)

# F-distribution
curve(df(x, df1 = 1, df2 = 1), from = 0, to = 4, n = 5000, col= col.b, lwd=lwd, ylim = c(0,1), xlab = "", ylab = "Density", main = "F-distribution")
curve(df(x, df1 = 10, df2 = 20), add = TRUE, n = 5000, col= col, lwd=lwd)
curve(df(x, df1 = 20, df2 = 1), add = TRUE, n = 5000, col= col.r, lwd=lwd)
abline(h=seq(0,1, by = .1), lty = 3, lwd = .3)
title(xlab = "x",line=2.2, cex.lab=1.2) 
legend("topright",legend = c("df1 = 1, df2 = 1", "df1 = 10, df2 = 20", "df1 = 20, df2 = 1"),lty = c(1,1,1), col =c(col.b, col, col.r), lwd = 2)

# Normal1
curve(expr = dnorm(x = x, mean=0,sd=1), from = -5, to = 5, col=col.b, lwd=lwd, ylim = c(0,1), xlab ="", ylab = "Density", main = "Normal")
curve(expr = dnorm(x = x, mean=0,sd=2), from = -5, to = 5, col=col, lwd=lwd, ylim = c(0,1), add = T)
curve(expr = dnorm(x = x, mean=2,sd=1), from = -5, to = 5, col=col.r, lwd=lwd, ylim = c(0,1), add = T)
abline(h=seq(0,1, by = .1), lty = 3, lwd = .3)
title(xlab = "x",line=2.2, cex.lab=1.2) 
legend("topright",legend = c("Normal, m = 0, sd = 1","Normal, m = 0, sd = 2","Normal, m = 2, sd = 1"),lty = c(1,1,1), col =c(col.b,col,col.r), lwd = 2)

# Normal2
curve(dt(x, df=1), from=-5, to=5, col=col.b, lty = 1, lwd=lwd, ylim = c(0,1), xlab ="", ylab = "Density", main = "t-distribution")
curve(expr = dnorm(x = x, mean=0,sd=1), from = -5, to = 5, col=col, lwd=lwd, ylim = c(0,1), add=TRUE)
curve(dt(x, df=10), from=-5, to=5, col=col.r, lty = 2, lwd=lwd, ylim = c(0,1),add = TRUE)
abline(h=seq(0,1, by = .1), lty = 3, lwd = .3)
title(xlab = "x",line=2.2, cex.lab=1.2) 
legend("topright",legend = c("t-distribution, df = 1", "t-distribution, df = 10", "Normal, m = 0, sd = 1"),lty = c(1,2,1), col =c(col.b,col.r,col), lwd = 2)

# Log normal 
curve(dlnorm(x, meanlog=0, sdlog=1), from=0, to=10, col=col.b, lwd=lwd, ylim = c(0,1), xlab = "", ylab = "Density", main = "Log-normal")
curve(dlnorm(x, meanlog=1.5, sdlog=1.2), add = TRUE, col=col, lwd=lwd, ylim = c(0,1))
curve(dlnorm(x, meanlog=1.8, sdlog=.15), add = TRUE, col=col.r, lwd=lwd, ylim = c(0,1))
abline(h=seq(0,1, by = .1), lty = 3, lwd = .3)
title(xlab = "x",line=2.2, cex.lab=1.2) 
legend("topright",legend = c("meanlog = 0, sdlog = 1", "meanlog = 1.5, sdlog = 1.2", "meanlog = 1.8, sdlog = 0.15"),lty = c(1), col =c(col.b, col, col.r), lwd = 2)

# Logistic
curve(dlogis(x,location = 0, scale = 1), from=-10, to=10, col = col.b, lwd=lwd, xlab ="",ylim = c(0,1), ylab = "Density", main = "Logistic")
curve(dlogis(x,location = 2, scale = 2), add = TRUE, col = col, lwd=lwd)
abline(h=seq(0,1, by = .1), lty = 3, lwd = .3)
title(xlab = "x",line=2.2, cex.lab=1.2) 
legend("topright",legend = c("loc = 0, scale = 1", "loc = 2, scale = 2"),lty = c(1), col =c(col.b, col), lwd = 2)

# unifrom
curve(dunif(x, min = 8,max = 9), from=5, to=25, col=col.b, lwd=lwd, ylim = c(0,1), xlab ="",ylab = "Density", main = "Uniform continuous") 
curve(dunif(x, min = 10,max = 15), from=5, to=25, col=col, lwd=lwd, ylim = c(0,1), add=T) 
curve(dunif(x, min = 16,max = 18), from=5, to=25, col=col.r, lwd=lwd, ylim = c(0,1), add = T) 
abline(h=seq(0,1, by = .1), lty = 3, lwd = .3)
title(xlab = "x",line=2.2, cex.lab=1.2) 
legend("topright",legend = c("min = 8, max = 9","min = 10, max = 15","min = 16, max = 18"),lty = c(1,1,1), col =c(col.b,col,col.r), lwd = 2)


## ----similar_plot_from_different_distributions, echo=FALSE, fig.width=7,fig.height=3--------------------------------------------------------------------------
par(mfrow = c(1,2))
set.seed(1234)
hist(rbinom(10000, 10, 0.5), #breaks = seq(-0.5, 10.5, by = 1), 
     main = "Binomial")
hist(rnorm(10000, 5, 1.5), xlim = c(0,10), main = "Normal")


## ----equivalence_between_distributions_Bern_binom, echo=FALSE, fig.width=7,fig.height=4-----------------------------------------------------------------------
par(mfrow = c(1,2), cex = 1.4)
# install.packages("Rlab")
library("Rlab")
# Bernoulli 
x <- seq(0, 1, by = 1) 
y_dbern <- dbern(x, prob = 0.7)
barplot(height = y_dbern,
        names.arg = setNames(c(.3, .7), c('absent (0)', 'present (1)')), 
        ylim = c(0, 1), xlab = '', ylab = 'probability', main = 'Bernoulli p = 0.7')

# Binomial (to get Bernoulli) 
x <- seq(0, 1, by = 1) 
y_binom <- dbinom(x,size = 1, prob = 0.7)
barplot(height = y_binom,
        names.arg = setNames(c(.3, .7), c('absent (0)', 'present (1)')), 
        ylim = c(0, 1), xlab = '', ylab = 'probability', main = 'Binomial, n = 1, p = 0.7')

# Simulate 10 (fair) coin flips
# set.seed(98765)
# rbinom(n = 10, size = 1, prob = 0.5)


## ----simulate_coin_flips_plot, echo=FALSE, eval=FALSE---------------------------------------------------------------------------------------------------------
## # Coin flips visualize
## set.seed(1235)
## 
## nb.flips = 10
## set.seed(98765)
## coin.flips = rbinom(n = nb.flips, size = 1, prob = 0.5)
## heads.tails = ifelse(coin.flips==1,"H","T")
## radius =1
## # initialize a plot
## plot(x = c(0, nb.flips*(radius+.5)), y = c(-.5, 1), type = "n",
##      axes = F,
##      asp = 1,
##      ylab = "", xlab = "")
## w = 0
## pos=0
## for (i in 1:nb.flips) {
##   if (i %% 5 == 0) {
##     pos = pos + 2.5
##     w = 0
##   }
## 
##   # prepare "circle data"
##   radius = 1
##   center_x = w + 1
##   center_y = pos
##   theta = seq(0, 2 * pi, length = 200) # angles for drawing points around the circle
## 
##   # draw the circle
##   lines(x = radius * cos(theta) + center_x, y = radius * sin(theta) + center_y)
##   polygon(x = radius * cos(theta) + center_x, y = radius * sin(theta) + center_y, col = scales::alpha("beige",.5))
##   text(center_x,center_y,labels = heads.tails[i], cex = 3)
##   w = w + 2*radius+.5
## }
## 


## ----dbinmon_prob.5-------------------------------------------------------------------------------------------------------------------------------------------
dbinom(2,3,prob = .5) # 0.375


## ----dbinom_plot_ex, echo=-c(1:5), fig.width=8,fig.height=4---------------------------------------------------------------------------------------------------
par(mfrow = c(1,1), mar = c(4,4,1,1))
set.seed(12345); lwd = 3
col = scales::alpha("black",.5); col.r = scales::alpha("red",.5)
col.g = scales::alpha("green",.5); col.b = scales::alpha("blue",.8)
add.l <- function(by=.1) {  abline(h=seq(0,1, by = by), lty = 3, lwd = .3)}
n = 7; x = 6; prob = 0.5
plot(0:8, dbinom(x = 0:8, size=n, prob=0.5),type='h', col = col.b, lwd=lwd, ylim = c(0,.4), ylab = "Probability", main = "s=7 and s=5; p=0.5 "); add.l(.05)
points(x = x, y = dbinom(x = x, size=n, prob=0.5), pch =19, col = col.b)

lines(c(0:8)+0.1, dbinom(x = 0:8, size=5, prob=0.5),type='h', col = col.r, lwd=lwd)
points(x = 3+.1, y = dbinom(x = 3, size=5, prob=0.5), pch =19, col = col.r)


## ----tables_rbinom_ex-----------------------------------------------------------------------------------------------------------------------------------------
set.seed(12345)
# 1000 experiments where each time, e.g., flipping a coin, I can either have a
table(rbinom(n = 1000, size=1, prob=0.5), dnn=NULL) # success or failure with p=.5
# 1 experiment where I have 1000 coins Where I sum all successes with p=.5
table(rbinom(n = 1, size=1000, prob=0.5), dnn=NULL) 
# 1000 experiments where each time, for example flipping 10 coins, where I 
table(rbinom(n = 1000, size=10, prob=0.5), dnn=NULL) # sum the success with p=.5


## ----Statistical_dist_binom, echo=-c(1:6), fig.width=8,fig.height=4-------------------------------------------------------------------------------------------
par(mfrow = c(2,3), mar = c(4,4,1,1))
set.seed(12345)
col = scales::alpha("black",.5); col.r = scales::alpha("red",.5)
col.g = scales::alpha("green",.5); col.b = scales::alpha("blue",.8)
lwd=2
add.l <- function(by=.1) {  abline(h=seq(0,1, by = by), lty = 3, lwd = .3)}
# For dbinom, x is the vector of quantiles 
plot(0:1, dbinom(x = 0:1, size=1, prob=0.5),type='h', col = col.b, lwd=lwd, ylim = c(0,1), ylab = "Probability", main = "s=1, p=.5"); add.l()
plot(0:20, dbinom(x = 0:20, size=20, prob=0.4),type='h', col = col.b, lwd=lwd, ylim = c(0,.2), ylab = "Probability", main = "s=20, p=.4");add.l()
plot(0:23, dbinom(x = 0:23, size=20, prob=0.9),type='h', col = col.b, lwd=lwd, ylim = c(0,0.3), ylab = "Probability", main = "s=20, p=.9");add.l()
plot(40:110, dbinom(x = 40:110, size=150, prob=0.5),type='h', col = col.b, lwd=lwd, ylim = c(0,.15), ylab = "Probability", main = "s=150, p=.5");add.l()


# hist(rbinom(n = 100000, size=50, prob=0.5), xlab = "Nb of successes", breaks = 100, col = col.b, main = "s=50, p=.5", probability = F)
plot(40:51, dbinom(x = 40:51, size=50, prob=0.99),type='h', col = col.b, lwd=lwd, ylim = c(0,.8), ylab = "Probability", main = "s=50, p=.99");add.l()
hist(rbinom(n = 100000, size=50, prob=0.99),xlab = "Nb of successes", breaks = 100, col = col.b, main = "s=50, p=.99", probability = F)


## ----binom_chall_coin, echo=-1, fig.width=8,fig.height=4------------------------------------------------------------------------------------------------------
par(mar =c(4,4,1,1))
mybin <- function(x,n,p) { choose(n, x)* p^x* (1-p)^(n-x) } 
mybin(0,4,.5) + mybin(4,4,.5)
dbinom(x = 0,4,.5) + dbinom(x = 4,4,.5)
lwd = 5
curve(dbinom(x,4,.5),0,6, n = 7,type = "h", lwd = lwd)
curve(dbinom(x,4,.5),0,0, type = "h",add = TRUE, col = "red", lwd = lwd)
curve(dbinom(x,4,.5),4,4, type = "h",add = TRUE, col = "red", lwd = lwd)

0.5^4*2


## ----petri_dish_bacteria_colonny, echo=FALSE, fig.width=5,fig.height=3----------------------------------------------------------------------------------------
set.seed(1235)
# initialize a plot
plot(c(-1, 3.5), c(-1, 1), type = "n", axes = F, asp = 1, ylab = "", xlab = "")
n = 25
# prepare "circle data"
radius = 1
center_x = 0
center_y = 0
theta = seq(0, 2 * pi, length = 200) # angles for drawing points around the circle

# draw the circle
lines(x = radius * cos(theta) + center_x, y = radius * sin(theta) + center_y)
polygon(x = radius * cos(theta) + center_x, y = radius * sin(theta) + center_y, col = scales::alpha("beige",.5))

center_x = 2.5
center_y = 0
lines(x = radius * cos(theta) + center_x, y = radius * sin(theta) + center_y)
polygon(x = radius * cos(theta) + center_x, y = radius * sin(theta) + center_y, col = scales::alpha("beige",.5))

rdmunif<-runif(n,0,1)
r = runif(n,0,radius) #radius * sqrt(rdmunif)
theta = rdmunif * 2 * pi
x = center_x + sqrt(r*radius) * cos(theta)
y = center_y + sqrt(r*radius) * sin(theta)

points(x,y, pch =19, cex = .6)

arrows(1.1,0,1.4,0, length = .1)
# draw sq
x.pos = c(1.9,2.5, 2.9,2,2.2)
y.pos = c(0,.3,-.6,.3,-.8)
dx = .1
for (i in 1:length(x.pos)) {
  polygon(x = c(x.pos[i],x.pos[i]+dx,x.pos[i]+dx,x.pos[i]),y = c(y.pos[i],y.pos[i],y.pos[i]-dx,y.pos[i]-dx))
}


## ----rpois_example_1.2----------------------------------------------------------------------------------------------------------------------------------------
set.seed(5937); rpois(n = 1,lambda = 2)


## ----Statistical_dist_poisson_wolf, echo=FALSE, fig.width=8,fig.height=2--------------------------------------------------------------------------------------
par(mfrow = c(1,1), mar = c(4,4,1,1))
set.seed(12345)
col = scales::alpha("black",.5); col.r = scales::alpha("red",.5)
col.g = scales::alpha("green",.5); col.b = scales::alpha("blue",.8)
lwd=2
add.l <- function(by=.1) {  abline(h=seq(0,1, by = by), lty = 3, lwd = .3)}
plot(0:10, dpois(x = 0:10, lambda=2),type='h', col = col.b, lwd=lwd, ylim = c(0,.41), ylab = "Probability", main = "lambda = 2"); add.l()


## ----Statistical_dist_poisson, echo=-c(1:6), fig.width=8,fig.height=4-----------------------------------------------------------------------------------------
# Poisson
par(mfrow = c(2,3), mar = c(4,4,1,1))
set.seed(12345)
col = scales::alpha("black",.5); col.r = scales::alpha("red",.5)
col.g = scales::alpha("green",.5); col.b = scales::alpha("blue",.8)
lwd=2
add.l <- function(by=.1) {  abline(h=seq(0,1, by = by), lty = 3, lwd = .3)}
# For dbinom, x is the vector of quantiles 
plot(0:5, dpois(x = 0:5, lambda=0.05),type='h', col = col.b, lwd=lwd, ylim = c(0,1), ylab = "Probability", main = "lambda = 0.05");add.l();abline(v = 0.05, lty = 3)
plot(0:10, dpois(x = 0:10, lambda=1),type='h', col = col.b, lwd=lwd, ylim = c(0,1), ylab = "Probability", main = "lambda = 1"); add.l();abline(v = 1, lty = 3)
plot(0:20, dpois(x = 0:20, lambda=5),type='h', col = col.b, lwd=lwd, ylim = c(0,0.3), ylab = "Probability", main = "lambda = 5");add.l();abline(v = 5, lty = 3)
plot(0:40, dpois(x = 0:40, lambda=20),type='h', col = col.b, lwd=lwd, ylim = c(0,.2), ylab = "Probability", main = "lambda = 20");add.l();abline(v = 20, lty = 3)


hist(rpois(n = 1000, lambda=25), xlab = "Nb of successes", breaks = 100, col = col.b, main = "n = 1000, lambda = 25", probability = F);abline(v = 25, lty = 3)
hist(rpois(n = 100000, lambda=50),xlab = "Nb of successes", breaks = 100, col = col.b, main = "n = 100000, lambda = 50", probability = F);abline(v = 50, lty = 3)


## ----equivalence_between_distributions_Binom_poisson, echo=FALSE----------------------------------------------------------------------------------------------
par(mfrow = c(2,2), mar = c(4,4,1,1))
x <- 0:10
n <- 10000
barplot(dbinom(x, n, 2/n), names.arg = x, ylim = c(0, 0.35), main = paste("Binomial with n = ", n))
barplot(dbinom(x, n, 9/n), names.arg = x, ylim = c(0, 0.35), main = paste("Binomial with n = ", n))
barplot(dpois(x, 2), names.arg = x, ylim = c(0, 0.35), main = paste("Poisson with Lambda = ", 2))
barplot(dpois(x, 9), names.arg = x, ylim = c(0, 0.35), main = paste("Poisson with Lambda = ", 9))

pbinom(q = 2, size = n, prob = 2/n)
ppois(q = 2, lambda = 2)


## ----uniform_distcrete_fun_dist, echo=-1, fig.width=5,fig.height=3--------------------------------------------------------------------------------------------
par(mar = c(3,4,1,1))
# Define uniform discrete 
dunifdisc<-function(x, min=0, max=1) ifelse(x>=min & x<=max & round(x)==x, 1/(max-min+1), 0)
runifdisc<-function(n, min=0, max=1) sample(min:max, n, replace=T)
curve(dunifdisc(x, 7,10), type = "h",from=6, to=11, col="black", xlab = "",
      lwd=1, ylim = c(0,.5), ylab = "Probability", main = "Uniform discrete"); title(xlab = "x",line=2.2, cex.lab=1.2) 


## ----uniform_distcrete_rdm_nb, echo=TRUE----------------------------------------------------------------------------------------------------------------------
set.seed(12345)
sample(7:10, size = 1,replace = T)
runifdisc(1,7,10)


## ----uniform_cont_angles, echo=-1, fig.width=5,fig.height=3---------------------------------------------------------------------------------------------------
par(mar = c(3,4,1,1))
x = seq(0,360,by = .1)
y <- dunif(x, min = min(x),max = max(x))
plot(y~x, type = "n", ylim = c(0,1.5*max(y)), ylab = "Density", xlab = "", main = "Uniform continuous")
polygon(c(x, rev(x), 0), c(y, rep(0,length(y)), 0), col=scales::alpha("blue",.5))
title(xlab = "x",line=2.2, cex.lab=1.2)


## ----uniform_continuous_rdm_nb, echo=TRUE---------------------------------------------------------------------------------------------------------------------
set.seed(12345)
runif(n = 1, min = 0, max = 360)


## ----Iris_normal, echo=FALSE, fig.width=12,fig.height=4-------------------------------------------------------------------------------------------------------
par(mfrow = c(1,3), cex = 1.4)
spiris = unique(iris$Species)
for (i in 1:length(spiris)) {
  tmp.iris=iris[iris$Species %in% spiris[i],]
  titl = paste("Iris",spiris[i])
  hist(tmp.iris$Sepal.Length,
       main = bquote(italic(.(titl))),
       xlab = "Sepal length",
       xlim = c(3.8,9), breaks = 10)
}


## ----Density_normal, fig.width=5,fig.height=5-----------------------------------------------------------------------------------------------------------------
curve(expr = dnorm(x = x, mean=0,sd=1), 
      from = -5, to = 5, ylim = c(0,1),
      col="black", ylab = "Density",
      main = "Density normal") 
abline(h=seq(0,1, by = .1), 
       lty = 3, lwd = .3)



## ----Density_normal_dissection_all_info, echo=FALSE, fig.width=7,fig.height=6---------------------------------------------------------------------------------
par(mfrow = c(1,1), mar = c(4,4,2.5,1))
x <- seq(-5, 5, 0.1)
cex  = 1.0
plot(x, dnorm(x, 0, 1), 
     main = "Density normal", type = "l", lwd = 3, 
     col = "black", ylab = "", xlab = "x", ylim = c(0,.5))
# quantile.normal = qnorm(seq(0,1, by = .05))
quantile.normal = qnorm(c(.999,0.001,.975,.025, .95, .05, .9, .1))
abline(v = quantile.normal[!is.infinite(quantile.normal)], lty = 3, lwd = .3) 
abline(h=seq(0,1, by = .1), lty = 3, lwd = .3)

p = 0.025 

# add the polygon to the left  
lb <- min(x) # Lower bound
ub <- qnorm(p)   # Upper bound
x2 <- seq(min(x), ub, length = 100) # New Grid
y <- dnorm(x2, 0, 1) # Densitypolygon(c(lb, x2, ub), c(0, y, 0), col = rgb(0, 0, 1, alpha = 0.5))
polygon(c(lb, x2, ub), c(0, y, 0), col = rgb(0, 0, 1, alpha = 0.5))
text(x = ub, y = p+0.04,labels = paste0("pnorm(qnorm(p)) = \n",p*100,"%"),adj = 0,pos = 2, cex=cex)

text(x = ub-.6, y = .2,labels = paste0("qnorm(p) = ", round(qnorm(p),2)),adj = 0,pos = 1,offset = -1, cex=cex)
arrows(x0 = ub,x1 = ub, y0 = .2,y1 = .1,code = 2,length=.1)


# add the polygon to the right 
lb <- qnorm(1-p) # Lower bound
ub <- max(x)   # Upper bound
x2 <- seq(lb, max(x), length = 100) # New Grid
y <- dnorm(x2, 0, 1) # Density
polygon(c(lb, x2, ub), c(0,y,0), col = rgb(0, 0, 1, alpha = 0.5))
text(x = lb, y = p+0.04,labels = paste0(p*100,"%"),adj = 0,pos = 4, cex=cex)

text(x = lb+.7, y = .2,labels = paste0(round(qnorm(1-p),2)," = qnorm(1-p)"),adj = 0,pos = 1,offset = -1, cex=cex)
arrows(x0 = lb,x1 = lb, y0 = .2,y1 = .1,code = 2,length=.1)

# Add the middle (red) polygon 
lb <- qnorm(p) # Lower bound
ub <- qnorm(1-p)   # Upper bound
x2 <- seq(lb, ub, length = 100) # New Grid
y <- dnorm(x2, 0, 1) # Density
polygon(c(lb, x2, ub), c(0,y,0), col = rgb(1, 0, 0, alpha = 0.5))
text(x = mean(x2), y = .15,labels = paste0((1-2*p)*100,"% \n=", "\npnorm(qnorm(1-p)) - \npnorm(qnorm(p))"),adj = 0,pos = 1, cex=cex)

text(x = 0, y = .47,labels = paste0(quote("Mean")),adj = 0,pos = 1,offset = 0, cex=cex)
arrows(x0 = 0,x1 = 0, y0 = .45,y1 = .42,code = 2,length=.1)

legend("topright",legend = c("Density normal (dnorm)","Density rnorm(100)"), lwd = 1, lty = 1, col = c("black", "green"))

set.seed(123)
rndat =rnorm(100)
# mean(rndat)
dens.nor = density(rndat)
lines(dens.nor, lwd = 3,col = scales::alpha("green",.8))
mybins=hist(rndat, plot = F, density = T, breaks = 100)
crn = mybins$density
brn = mybins$breaks


## ----pnorm_qnorm2---------------------------------------------------------------------------------------------------------------------------------------------
pnorm(1.645) # Gives the probability at X 

qnorm(p = 0.05, lower.tail = F) # Gives X at a prob.

qnorm(p = 0.025, lower.tail = F) # Value for 5% in both tails


## ----normal_shade_95_5, echo=FALSE, fig.width=5,fig.height=5--------------------------------------------------------------------------------------------------
x <- seq(-5, 5, 0.1)
cex  = 1
plot(x, dnorm(x, 0, 1), main = "Density normal", type = "l", lwd = 3, col = "black", ylab = "", xlab = "x", ylim = c(0,1))
quantile.normal = qnorm(c(.95, .05))
abline(v = quantile.normal[!is.infinite(quantile.normal)], lty = 3, lwd = .3) 
abline(h=seq(0,1, by = .1), lty = 3, lwd = .3)

p = 0.05

# add the polygon to the right 
lb <- qnorm(1-p) # Lower bound
ub <- max(x)   # Upper bound
x2 <- seq(lb, max(x), length = 100) # New Grid
y <- dnorm(x2, 0, 1) # Density
polygon(c(lb, x2, ub), c(0,y,0), col = rgb(1, 0, 0, alpha = 0.5))
text(x = lb, y = p+0.04,labels = paste0(p*100,"%"),adj = 0,pos = 4, cex=cex)

# Add the middle (blue) polygon 
lb <- min(x) # Lower bound
ub <- qnorm(1-p)   # Upper bound
x2 <- seq(lb, ub, length = 100) # New Grid
y <- dnorm(x2, 0, 1) # Density
polygon(c(lb, x2, ub), c(0,y,0), col = rgb(0, 0, 1, alpha = 0.5))
text(x = 0, y = .2,labels = paste0((1-p)*100,"%"),adj = 0,pos = 1, cex=cex)

text(x = 0, y = .55,labels = paste0(quote("Mean")),adj = 0,pos = 1,offset = 0, cex=cex)
arrows(x0 = 0,x1 = 0, y0 = .5,y1 = .42,code = 2,length=.1)


## ----normal_area_function, echo=FALSE, eval=TRUE--------------------------------------------------------------------------------------------------------------
draw.normal <- function(mean = 0, sd = 1, set.seed=1, prob = 0.025, text = FALSE, text.height = .55, where = c("both","left","right","middle"), middle = c(-1,1)) {
  set.seed(set.seed)
  x <- seq(-5, 5, 0.1)
  cex = 1
  plot(x, dnorm(x, mean, sd), 
       # main = "Density normal", 
       main = "", 
       type = "l", lwd = 3, col = "black", 
       ylab = "", xlab = "",
       ylim = c(0,1))
  title(xlab = "x", ylab="", line=2.2, cex.lab=1.2)

  # abline(v = quantile.normal[!is.infinite(quantile.normal)], lty = 3, lwd = 1, col = c("black")) 
  
  # Horizontal 
  abline(h=seq(0,1, by = .1), lty = 3, lwd = .3)
  
  p = prob 
  
  if(where=="both"){
    # add the polygon to the left  
    lb <- min(x) # Lower bound
    ub <- qnorm(p)   # Upper bound
    x2 <- seq(min(x), ub, length = 100) # New Grid
    y <- dnorm(x2, 0, 1) # Densitypolygon(c(lb, x2, ub), c(0, y, 0), col = rgb(0, 0, 1, alpha = 0.5))
    polygon(c(lb, x2, ub), c(0, y, 0), col = rgb(0, 0, 1, alpha = 0.5))
    text(x = -2, y = .2,
         labels = paste0(p*100,"%"),adj = 0,pos = 2, cex=cex)
    
    if (text) {
      text(x = ub, y = .2,
           labels = paste0(round(qnorm(p),2)),adj = 0,pos = 1,offset = 0, cex=cex)
      arrows(x0 = ub,x1 = ub, y0 = .2,y1 = .1,code = 2,length=.1)
    }
    
    # add the polygon to the right 
    lb <- qnorm(1-p) # Lower bound
    ub <- max(x)   # Upper bound
    x2 <- seq(lb, max(x), length = 100) # New Grid
    y <- dnorm(x2, 0, 1) # Density
    polygon(c(lb, x2, ub), c(0,y,0), col = rgb(0, 0, 1, alpha = 0.5))
    text(x = 2, y = .2,
         labels = paste0(p*100,"%"),adj = 0,pos = 4, cex=cex)
    
    if (text) {
      text(x = lb, y = .2,
           labels = paste0(round(qnorm(1-p),2)),adj = 0,pos = 1,offset = 0, cex=cex)
      arrows(x0 = lb,x1 = lb, y0 = .2,y1 = .1,code = 2,length=.1)
    }
    
    # Add the middle (red) polygon 
    lb <- qnorm(p) # Lower bound
    ub <- qnorm(1-p)   # Upper bound
    x2 <- seq(lb, ub, length = 100) # New Grid
    y <- dnorm(x2, 0, 1) # Density
    polygon(c(lb, x2, ub), c(0,y,0), col = rgb(1, 0, 0, alpha = 0.5))
    text(x = mean(x2), y = text.height,
         labels = paste0((1-2*p)*100,"%"),adj = 0,pos = 1, cex=cex)
    
    if (text) {
      text(x = 0, y = .55,
           labels = paste0(quote("Mean")),adj = 0,pos = 1,offset = 0, cex=cex)
      arrows(x0 = 0,x1 = 0, y0 = .5,y1 = .42,code = 2,length=.1)
    }
  }
  
  if (where=="left") {
    # add the polygon to the left  
    lb <- min(x) # Lower bound
    ub <- qnorm(p)   # Upper bound
    x2 <- seq(min(x), ub, length = 100) # New Grid
    y <- dnorm(x2, 0, 1) # Densitypolygon(c(lb, x2, ub), c(0, y, 0), col = rgb(0, 0, 1, alpha = 0.5))
    polygon(c(lb, x2, ub), c(0, y, 0), col = rgb(0, 0, 1, alpha = 0.5))
    text(x = -2, y = .2,
         labels = paste0(p*100,"%"),adj = 0,pos = 2, cex=cex)
    
    if (text) {
      text(x = ub, y = .2,
           labels = paste0(round(qnorm(p),2)),adj = 0,pos = 1,offset = 0, cex=cex)
      arrows(x0 = ub,x1 = ub, y0 = .2,y1 = .1,code = 2,length=.1)
    }
    
    # Add the middle (red) polygon 
    lb <- qnorm(p) # Lower bound
    ub <- max(x)   # Upper bound
    x2 <- seq(lb, ub, length = 100) # New Grid
    y <- dnorm(x2, 0, 1) # Density
    polygon(c(lb, x2, ub), c(0,y,0), col = rgb(1, 0, 0, alpha = 0.5))
    text(x = 2, y = .2,
         labels = paste0((1-p)*100,"%"),adj = 0,pos = 4, cex=cex)
    
    if (text) {
      text(x = 0, y = .55,
           labels = paste0(quote("Mean")),adj = 0,pos = 1,offset = 0, cex=cex)
      arrows(x0 = 0,x1 = 0, y0 = .5,y1 = .42,code = 2,length=.1)
    } 
    
  }
  
  if (where=="right"){
    # add the polygon to the left  
    lb <- min(x) # Lower bound
    ub <- qnorm(1-p)   # Upper bound
    x2 <- seq(min(x), ub, length = 100) # New Grid
    y <- dnorm(x2, 0, 1) # Densitypolygon(c(lb, x2, ub), c(0, y, 0), col = rgb(0, 0, 1, alpha = 0.5))
    polygon(c(lb, x2, ub), c(0, y, 0), col = rgb(1, 0, 0, alpha = 0.5))
    text(x = -2, y = .2,
         labels = paste0((1-p)*100,"%"),adj = 0,pos = 2, cex=cex)
    
    if (text) {
      text(x = ub, y = .2,
           labels = paste0(round(qnorm(p),2)),adj = 0,pos = 1,offset = 0, cex=cex)
      arrows(x0 = ub,x1 = ub, y0 = .2,y1 = .1,code = 2,length=.1)
    }
    
    # Add the middle (red) polygon 
    lb <- qnorm(1-p) # Lower bound
    ub <- max(x)   # Upper bound
    x2 <- seq(lb, ub, length = 100) # New Grid
    y <- dnorm(x2, 0, 1) # Density
    polygon(c(lb, x2, ub), c(0,y,0), col = rgb(0, 0, 1, alpha = 0.5))
    text(x = 2, y = .2,
         labels = paste0((p)*100,"%"),adj = 0,pos = 4, cex=cex)
    
    if (text) {
      text(x = 0, y = .55,
           labels = paste0(quote("Mean")),adj = 0,pos = 1,offset = 0, cex=cex)
      arrows(x0 = 0,x1 = 0, y0 = .5,y1 = .42,code = 2,length=.1)
    }
    
  }
  
    if (where=="middle"){
    # Add the middle (red) polygon 
    lb <- (middle[1]) # Lower bound
    ub <- (middle[2])   # Upper bound
    x2 <- seq(lb, ub, length = 100) # New Grid
    y <- dnorm(x2, 0, 1) # Density
    polygon(c(lb, x2, ub), c(0,y,0), col = rgb(0, 0, 1, alpha = 0.5))
    text(x = mean(middle), y = dnorm(mean(middle))+.2,
         labels = paste0(round(pnorm(middle[2])-pnorm(middle[1]),digits = 4)*100,"%"),adj = 0,pos = 1, cex=cex)
    
  }

}


## ----normal_dist_area, echo=-1, fig.width=9,fig.height=5------------------------------------------------------------------------------------------------------
par(mfrow=c(2,3), mar = c(4,4,.1,.1), cex = 1.1)
# The function is not shown, but can be found in the script (Markdown)
draw.normal(where = "both",  prob = 0.05/2)
draw.normal(where = "both",  prob = 0.2/2 )
draw.normal(where = "both",  prob = 0.5/2 )
draw.normal(where = "both",  prob = 0.95/2)
draw.normal(where = "left",  prob = 0.05  )
draw.normal(where = "right", prob = 0.05  )


## ----Normal_pdf_important_values------------------------------------------------------------------------------------------------------------------------------
sd = 1
probability.left.side = (pnorm(q = c(sd*1,sd*2,sd*3),lower.tail = F)*100)
probability.right.side = (pnorm(q = c(sd*1,sd*2,sd*3),lower.tail = T)*100)
percent.data.under.curve = probability.right.side - probability.left.side
p.from.mean = round(percent.data.under.curve,2)


## ----normal_dist_area2, echo=FALSE, fig.width=11,fig.height=3-------------------------------------------------------------------------------------------------
par(mfrow=c(1,3), mar = c(4,4,1,1), cex = 1.1)
draw.normal(where = "both",  prob = round(probability.left.side[1]/100,3))
draw.normal(where = "both",  prob = round(probability.left.side[2]/100,3))
draw.normal(where = "both",  prob = round(probability.left.side[3]/100,3))


## ----Normal_pdf_important_values2-----------------------------------------------------------------------------------------------------------------------------
qnorm(p = c(.75, .95,.975, .995), mean = 0, sd = 1, lower.tail = F) # notice the "lower.tail"
qnorm(p = c(.75, .95,.975, .995), mean = 0, sd = 1, lower.tail = T)


## ----normal_dist_area3, echo=FALSE, fig.width=14,fig.height=3-------------------------------------------------------------------------------------------------
par(mfrow=c(1,4), mar = c(4,4,1,1), cex = 1.1)
length.arrow.head = 0.10
draw.normal(where = "both",  prob = .25); x.val = qnorm(p = c(.25),lower.tail = T)
arrows(x0 = x.val-0.5, x1 = x.val, y0 = .5,y1 = .35, length = length.arrow.head)
text(-1.4, y = .6, labels = round(x.val,2))

draw.normal(where = "both",  prob = .05); x.val = qnorm(p = c(.05),lower.tail = T)
arrows(x0 = x.val, x1 = x.val, y0 = .5,y1 = .2, length = length.arrow.head)
text(x.val, y = .6, labels = round(x.val,2))

draw.normal(where = "both",  prob = .025); x.val = qnorm(p = c(.025),lower.tail = T)
arrows(x0 = x.val, x1 = x.val, y0 = .5,y1 = .15, length = length.arrow.head)
text(x.val, y = .6, labels = round(x.val,2))

draw.normal(where = "both",  prob = .005); x.val = qnorm(p = c(.005),lower.tail = T)
arrows(x0 = x.val+1, x1 = x.val, y0 = .5,y1 = .07, length = length.arrow.head)
text(x.val+1, y = .6, labels = round(x.val,2))



## ----normal_statistics_areas_1, echo=FALSE, fig.width=7,fig.height=5------------------------------------------------------------------------------------------
par(mfrow=c(1,1), mar = c(4,4,0,1))
hypothesis.testing <- function(mean.pop=2, mean.sn = 0,prob = .05, sd.sn =1, sd = 1, n=NULL) {
  
if(!is.null(n)){
  se = sd/sqrt(n)
  sd = se
  # print(se)
}
  
set.seed=1;cex =1
mean = mean.sn
p = prob
x = -6:10
mean2 = mean.pop
col1=rgb(0, 0, 1, alpha = 0.5)
col2=rgb(0, 1, 0, alpha = 0.5)
col3=rgb(0.8, 0.4, 0.2, alpha = 0.5)

dval=dnorm(x, mean = mean, sd = sd.sn)
dval2= dnorm(x, mean = mean2, sd = sd)
curve(dnorm(x, mean = mean, sd = sd.sn),-6,10, n = 10000,xlim = c(-5, 7), 
      ylim = range(0,1,max(dval,dval2)),
      xlab = "", lwd =3 )
# dnorm can have values that are GREATER THAN 1 # see https://stackoverflow.com/questions/42661973/r-density-plot-y-axis-larger-than-1
text(x = -4, y = .25,labels = "Sampling \ndistribution \nif H0 is true",adj = 0,pos = 4, cex=cex)
arrows(x0 = -3,x1 = -2, y0 = .2,y1 = .1,code = 2,length=.1)

null.test =qnorm(p = p, mean = mean, sd=sd.sn,lower.tail = F)
abline(v = null.test)



# Add the right (solid) polygon 
lb <- qnorm(pnorm(null.test,mean = mean,sd = sd.sn),mean,sd = sd.sn) # Lower bound
ub <- max(x)   # Upper bound
x2 <- seq(lb, ub, length = 100) # New Grid
y <- dnorm(x2, mean = mean, sd = sd.sn) # Density
polygon(c(lb, x2, ub), c(0,y,0), col = col3,
        # density = 10, angle = -45, 
        lwd = 2)

curve(dnorm(x, mean = mean2, sd = sd),-6,10, n = 10000, xlim = c(-5, 8), ylim = c(0,0.6),col = "red",add = T, lwd =3 )
text(x = 4, y = .25,labels = "Sampling \ndistribution \nof population",adj = 0,pos = 4, cex=cex)
arrows(x0 = 5,x1 = 4, y0 = .2,y1 = .1,code = 2,length=.1)

# add the polygon to the left  
lb <- min(x) # Lower bound
ub <- qnorm(pnorm(null.test,mean = mean2, sd=sd),mean = mean2,sd=sd)   # Upper bound
x2 <- seq(min(x), ub, length = 100) # New Grid
y <- dnorm(x2, mean2, sd) # Densitypolygon(c(lb, x2, ub), c(0, y, 0), col = rgb(0, 0, 1, alpha = 0.5))
polygon(c(lb, x2, ub), c(0, y, 0), col = col1,density = 10, angle = 45, lwd = 2)
# text(x = -2, y = .2,labels = paste0(p*100,"%"),adj = 0,pos = 2, cex=cex)

# Add the right (green) polygon 
lb <- qnorm(pnorm(null.test,mean = mean2, sd=sd),mean = mean2,sd=sd) # Lower bound
ub <- max(x)   # Upper bound
x2 <- seq(lb, ub, length = 100) # New Grid
y <- dnorm(x2, mean2, sd) # Density
polygon(c(lb, x2, ub), c(0,y,0), col = col2,density = 16, angle = -45, lwd = 2)
# text(x = 2, y = .2,labels = paste0((1-p)*100,"%"),adj = 0,pos = 4, cex=cex)


legend("topleft",
       legend = c("Type 1 error", "Type 2 error", "Power"),
       density = c(NA,10,16),
       angle = c( NA, 45, -45), 
       col =c(col3, col1,col2),
       fill = c(col3,col1,col2),
       ncol = 1,
       cex = 1, bg = "white"
)
# text(x = 2, y = .2,labels = paste0((1-p)*100,"%"),adj = 0,pos = 4, cex=cex)

}
hypothesis.testing(mean.pop = 1)


## ----normal_statistics_areas_2, echo=FALSE, fig.width=7,fig.height=5------------------------------------------------------------------------------------------
par(mfrow=c(1,1), mar = c(4,4,0,1), cex = 1.1)
hypothesis.testing(mean.pop = 2)


## ----normal_statistics_areas_3, echo=FALSE--------------------------------------------------------------------------------------------------------------------
par(mfrow=c(1,1), mar = c(4,4,0,1), cex = 1.1)
hypothesis.testing(2,sd = 5, n = 10)


## ----normal_statistics_areas_4, echo=FALSE--------------------------------------------------------------------------------------------------------------------
par(mfrow=c(1,1), mar = c(4,4,0,1), cex = 1.1)
hypothesis.testing(2,sd = 5, n = 100)


## ----Standard_normal_transformation, echo=FALSE, fig.width=10,fig.height=5------------------------------------------------------------------------------------
par(mfrow = c(1,2))
set.seed(1235)
x <- seq(-5, 20, 0.1)
cex  = 1
black.5 = scales::alpha("black",.5)
red.5 = scales::alpha("red",.5)
blue.5 = scales::alpha("blue",.5)
# Add the standard normal 
plot(x, dnorm(x, 0, 1), main = "Density normal", type = "l", lwd = 3, col = black.5, ylab = "", xlab = "x", ylim = c(0,1))
text(x = 0, y = .62,labels = paste0("Standard \nnormal"),adj = 0,pos = 1,offset = 0, cex=cex)
arrows(x0 = 0,x1 = 0, y0 = .5,y1 = .42,code = 2,length=.1)

# Add a population distribution 
lines(x, dnorm(x, 15, 2), main = "Density normal", type = "l", lwd = 3, col = red.5, ylab = "", xlab = "x", ylim = c(0,1))
text(x = 15, y = .55,labels = paste0("Normal \n mean=15 \n sd=2"),adj = 0,pos = 1,offset = 0, cex=cex)
arrows(x0 = 15,x1 = 15, y0 = .35,y1 = .25,code = 2,length=.1)

# Add simulated data 
my.data = rnorm(100, 15, 2)
mean.data = mean(my.data)
sd.data = sd(my.data)
lines(density(my.data), col = blue.5, lwd = 3)
# Add guides
abline(h=seq(0,1, by = .1), lty = 3, lwd = .3)

# Add the standard normal 
plot(x, dnorm(x, 0, 1), main = "Density normal", type = "l", lwd = 3, col = black.5, ylab = "", xlab = "x", ylim = c(0,1))

# Add a population distribution 
lines(x, dnorm(x, 0, 1), main = "Density normal", type = "l", lwd = 3, col = red.5, ylab = "", xlab = "x", ylim = c(0,1))
lines(density((my.data-mean.data)/sd.data), col = blue.5, lwd = 3)

# Add guides
abline(h=seq(0,1, by = .1), lty = 3, lwd = .3)



## ----rnom_function, eval=FALSE--------------------------------------------------------------------------------------------------------------------------------
## rnorm()


## ----rnom_function_example------------------------------------------------------------------------------------------------------------------------------------
set.seed(1234)
n <-10
rnorm(n)


## ----risland_function_example, fig.width=7,fig.height=5-------------------------------------------------------------------------------------------------------
size = 1000000
prob = 1/1000000
x = seq(0,10,by = 1)
plot(dbinom(x = x, size = size, prob = prob)~x, 
     type = "h", lwd = 4, ylab = "Probability", xlab = "Event")
points(dbinom(x = x[-1], size = size, prob = prob)~x[-1], 
       type = "h", lwd = 4, col = "red")

# This is the probability 
pbinom(0,size = size,prob = prob,lower.tail = F)
# pbinom(0,size = size,prob = prob,lower.tail = T)


## ----clt_challenger, eval=FALSE-------------------------------------------------------------------------------------------------------------------------------
## set.seed(42)
## sample(1:6, size = 1, replace = TRUE)


## ----clt_challenger_ans, echo=c(-1), fig.width=11,fig.height=3.5----------------------------------------------------------------------------------------------
par(mfrow = c(1,4),mar = c(4,4,1,1), cex = 1.2)
set.seed(12345)
roll = 1e4 # Nb of times we do the experiment (characterize the underlying r.v.)
for (nb.dice in c(1,2,10,50)) {
  res.exp = replicate(roll, simplify = T,
                      sample(1:6, size = nb.dice, replace = TRUE))
  if(nb.dice == 1){Xsum = res.exp } 
  else {Xsum = apply(res.exp, 2, sum)} # add independent r.v. (sum of cols)
  hist(Xsum, main=paste("n =",roll,ifelse(nb.dice==1,"dice =","die ="),nb.dice))}


## ----clt_challenger_ans_normalized, echo=c(-1), fig.width=9,fig.height=3.8------------------------------------------------------------------------------------
par(mfrow = c(1,4),mar = c(4,4,1,1), cex = 1.2)
set.seed(12345)
roll = 1e4 # Nb of times we do the experiment (characterize the underlying r.v.)
for (nb.dice in c(1,2,10,500)) {
  res.exp = replicate(roll, simplify = T,
                      sample(1:6, size = nb.dice, replace = TRUE))
  if(nb.dice == 1){Xsum = res.exp; Xmean = (res.exp)} else {
    Xsum = apply(res.exp, 2, sum) # add independent r.v. (sum of cols)
    Xmean = apply(res.exp, 2, mean)}
  hist((Xmean-mean(Xmean))/(sd(Xmean)), xlim = c(-4,4), ylim = c(0,1),
       main=paste("n =",roll,ifelse(nb.dice==1,"dice =","die ="),nb.dice), probability = T)
curve(dnorm, from = -10, to =10, add =TRUE)}


## ----simulation_clt_1, echo=-1, fig.width=8,fig.height=4.5----------------------------------------------------------------------------------------------------
par(mfrow=c(2,2), cex =1.1, mar = c(4,4,1,1)) # set window 
n = 1000 # Number of points 
# Generate multiple additions of random variables 
for(i in c(2, 50, 1000, 5000)){
  clt = replicate(i, rexp(n, rate = 1), simplify = FALSE)
  hist(apply(do.call(cbind,clt),1,sum), main = paste("Hist. of",i,"variables"), xlab = "x") # Draw the histogram 
}


## ----uniform_distcrete_fun, echo=FALSE------------------------------------------------------------------------------------------------------------------------
# Define uniform discrete 
dunifdisc<-function(x, min=0, max=1) ifelse(x>=min & x<=max & round(x)==x, 1/(max-min+1), 0)
punifdisc<-function(q, min=0, max=1) ifelse(q<min, 0, ifelse(q>=max, 1, (floor(q)-min+1)/(max-min+1)))
qunifdisc<-function(p, min=0, max=1) floor(p*(max-min+1))
runifdisc<-function(n, min=0, max=1) sample(min:max, n, replace=T)


## ----CLT_uniform_example_dice, echo=-1, eval=FALSE, fig.width=8,fig.height=5----------------------------------------------------------------------------------
## par(mfrow=c(2,2)) # set window
## # Generate multiple additions of random variables
## for(i in c(2, 50, 1000, 5000)){
##   clt = replicate(i, runifdisc(n),simplify = FALSE)
##   hist(apply(do.call(cbind,clt),1,sum), main = paste("Histogram of",i,"variables"),xlab = "x") # Draw the histogram
## }


## ----CLT_uniform_example_dice_hist, echo=FALSE, fig.width=8,fig.height=5--------------------------------------------------------------------------------------
par(mfrow=c(2,2)) # set window 
# Generate multiple additions of random variables 
for(i in c(1, 50, 1000, 5000)){
  clt = replicate(i, runifdisc(n,1,6),simplify = FALSE)
  hist(apply(do.call(cbind,clt),1,sum), main = paste("Histogram of",i,"variables"),xlab = "x") # Draw the histogram 
}


## ----uniform_disc_example_plot, echo=FALSE, fig.width=3,fig.height=3------------------------------------------------------------------------------------------
curve(dunifdisc(x, 7,10), type = "h",from=6, to=11, col="black", lwd=1, ylim = c(0,1), ylab = "Density", main = "Uniform") 


## ----clt_unif_discrete_plots, echo=-1, fig.width=8,fig.height=3-----------------------------------------------------------------------------------------------
par(mfrow=c(1,2), cex = 1.1) # set window 
n = 1e6 # Number of points 
# Generate multiple additions of random variables 
clt = replicate(n = 20, # flip 20 dice
                runifdisc(n = n, min = 1, max = 6), simplify = FALSE)
sum.rdm.var = apply(do.call(what = cbind, args = clt), 1, sum) 
mean.rdmv = mean(sum.rdm.var); sd.rdmv = sd(sum.rdm.var) # mean and sd
hist(sum.rdm.var, main = paste("Histogram of",20,"variables"), xlab = "x")
curve(expr = dnorm(x,mean = mean.rdmv,sd = sd.rdmv),
      from = mean.rdmv-5*sd.rdmv, to = mean.rdmv+5*sd.rdmv,ylab = "Density")


## ----clt_uniform_continuous.plots, echo=FALSE, fig.width=8,fig.height=4---------------------------------------------------------------------------------------
# modified from [Using R to simulate the Central Limit Theorem](https://consultglp.com/wp-content/uploads/2016/10/using-r-to-simulate-the-central-limit-theorem.pdf)
set.seed(2345)
# Simulation of central limit theorem #=================================================
layout(matrix(c(1,2,3,4),2,2,byrow=TRUE))
clt.fun <- function(size = 1, repeats = 10000,min = 0,max = 1, get.mean = TRUE) {
  v=runif(n = size*repeats, min = min,max = max) # Vector of uniform random variables. 
  w=matrix(data = v,nrow = size, ncol = repeats) # Enter v into a matrixsizeXrepeats). 
  dim(w)
  # w[1:4,1:4]
  y=colSums(w) # Sum the columns.
  
  if (get.mean) {
  y.mean=colSums(w)/size # Sum the columns.
  hist(y.mean,freq=FALSE,ann=FALSE, xlim = c(min, max)) # Histogram. 
  } else {
    hist(y,freq=FALSE,ann=FALSE) # Histogram.
    }
  title(paste("size",size))
  
}
#Sum of 1 uniform random variables simulated 10000 times 
clt.fun(size = 1)
#Sum of 2 uniform random variables simulated 10000 times 
clt.fun(size = 2)
#Sum of 4 uniform random variables simulated 10000 times 
clt.fun(size = 4)
#Sum of 20 uniform random variables simulated 10000 times 
clt.fun(size = 20)



## ----Galton_board, echo=FALSE---------------------------------------------------------------------------------------------------------------------------------
source(file = "scripts/Galton_board.R")
galtonqunincunx


## ----multivar_norm_pca, echo=-1, fig.width=5,fig.height=2-----------------------------------------------------------------------------------------------------
par(mfrow = c(1,2), mar = c(4,4,0.5,0.5)); b.5 = scales::alpha("black",0.5)
set.seed(24601) # setting this so the random results will be repeatable 
library(MASS)
covmat <- matrix(c(1.0,   0.2,   0.6, # variance covariance matrix of the data 
                   0.2,   2.0,  -0.5, # Positive-definite symmetric matrix
                   0.6,  -0.5,   1.0), nrow=3) 
data <- mvrnorm(n = 300,  
                mu = c(1,-1,0), # mean of the data 
                Sigma=covmat) # generate random data that match that variance covariance matrix
plot(data[,1:2], pch = 19, col = b.5); abline(h=0,v=0,lty = 3)
plot(data[,2:3], pch = 19, col = b.5); abline(h=0,v=0,lty = 3)


## ----multivar_norm_pca_plot, echo=-1, fig.width=6,fig.height=5------------------------------------------------------------------------------------------------
# biplot(prcomp(data)) # check var(data) to verify the covmat
ggplot2::autoplot(vegan::rda(data)) + theme_classic() + 
  theme(plot.margin = unit(c(0,0.5,0,0.5), "cm"), legend.position="none")+
  geom_hline(yintercept=0, linetype="dotted",color = "grey50", size=0.5)+
  geom_vline(xintercept=0, linetype="dotted",color = "grey50", size=0.5)



## ----multivar_norm_cov_and_correlation, echo=-1, fig.width=5,fig.height=2.5-----------------------------------------------------------------------------------
par(mfrow = c(1,1), mar = c(4,4,1.5,0.5)); b.5 = scales::alpha("black",0.5); xlim = c(8,14); ylim = c(7,14)
library(faux); set.seed(12345); n = 500 #number of values to generate for the variable x
traits = rnorm_multi(n = n, # ??rnorm_multi
                     mu = c(11.317531, 10.470744, 9.377967), # Mean of 3 variables
                     sd = c( 0.6247866, 0.7134755, 0.5502800), # SD of variables 
                     r = c(0.4887644, 0.4678446, 0.7056161), # correlation between variables 
                     varnames = c("Beak length", "Beak depth", "Beak width"), empirical = TRUE)
plot(traits[,1:2],col=b.5,pch=19,main="Simulated data",xlab="");title(xlab="Beak length",ylab="",line=2.2,cex.lab=1.2)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
cov(traits); cor(traits) # you can verify that these were the input in the function 


## ----simulation_PCA, echo=-1, fig.width=5,fig.height=2--------------------------------------------------------------------------------------------------------
par(mfrow = c(1,2), mar = c(4,4,0.5,0.5))
set.seed(123) # setting this so the random results will be repeatable 
library(MASS)
# Simulating 3 traits for 4 different species 
n = 200
Amat1 = MASS::mvrnorm(n, mu = c(11.2,11.8,9.91), Sigma = diag(c(1.31,1.01,1.02))) 
Amat2 = MASS::mvrnorm(n, mu = c(7.16,8.54,6.82), Sigma = diag(c(0.445,0.546,0.350)))
Amat3 = MASS::mvrnorm(n, mu = c(15.6,14.6,13.5), Sigma = diag(c(1.43,0.885,0.990)))
Amat4 = MASS::mvrnorm(n, mu = c(8.65,14.1,8.24), Sigma = diag(c(0.535,0.844,0.426)))
Amat = rbind(Amat1,Amat2,Amat3,Amat4)
Amat.gr = cbind(Amat, gl(4,k=n,labels = c(1,2,3,4)))

# by(Amat.gr[,1:3],INDICES = Amat.gr[,4],FUN = cov) # calc. cov. mat. for all gr 

summary(m1 <- prcomp(Amat, scale= T))


## ----simulation_PCA_plot, echo=-1, fig.width=5,fig.height=2---------------------------------------------------------------------------------------------------
par(mfrow = c(1,2), mar = c(4,4,0.5,0.5))
# biplot(m1, xlabs=rep(".", nrow(Amat)), cex = 3)
plot(Amat[,1],Amat[,2], pch = 19, col = gl(4,k=n,labels = c(1,2,3,4))); abline(h=mean(Amat[,1]),v=mean(Amat[,2]),lty = 3)
plot(vegan::scores(m1), asp = 1, pch = 19, col = gl(4,k=n,labels = c(1,2,3,4))); abline(h=0,v=0,lty = 3)
# library(ggvegan)
# autoplot(vegan::rda(Amat))


## ----parameter_distribution-----------------------------------------------------------------------------------------------------------------------------------
x = rnorm(100)
# x = rt(100,1000)
library(fitdistrplus)
descdist(x, discrete = FALSE)
fit.pois <- fitdist(x, "norm")
plot(fit.pois)


## ----lm_assumtion_expected_mean, echo=FALSE, fig.width=4,fig.height=4-----------------------------------------------------------------------------------------
par(mar=c(4,4,0.1,0.1))
set.seed(12345678)
n = 100; beta0 = 2.5; beta1 = 0.8
x.lm = rnorm(n = n, mean = 10, sd = 1)
err = rnorm(n = n, mean = 0, sd = 1)
# Linear combination 
y.lm = beta0 + beta1*x.lm + err
# Make a dataframe of the data 
df.lm = data.frame(x = x.lm, y = y.lm)
par(mar = c(4,4,.5,.5))
# Colour 
b.5 = scales::alpha("black",alpha = .5)

# PLot the data 
plot(y~x, data = df.lm, pch = 19, col = b.5, xlab = "")
title(xlab = "x", ylab="", line=2.2, cex.lab=1.2)
# Model the data 
lm.out = lm(y~x, data = df.lm)
# Add a line to the plot 
abline(lm.out)

####
slm=summary(lm.out)
# sd(slm$residuals)
yvals = seq(floor(range(df.lm$y)[1]),ceiling(range(df.lm$y)[2]),length.out = 1000)
xvals = seq(8,12, length.out = 10)
yh = dnorm(x =yvals, mean = predict(lm.out,newdata = data.frame(x = xvals)),sd = (slm$sigma))
# length(yh)
xrp = rep(xvals, length(yh)/length(xvals)) # each=length(yh)/length(xvals)
# length(xrp)
points(xrp,yvals, cex = yh/max(yh), pch = 19, col = scales::alpha("red",.5))


## ----lm_assumtion_expected_mean_interactive, echo=FALSE, fig.width=3,fig.height=3-----------------------------------------------------------------------------
library(plotly)
library(dplyr)

beta <- lm.out$coefficients
fxy <- function(x,y){  beta[1] + beta[2]*x }
# fxy2 <- function(x,y){  
#   observeddata <- cbind(predict(lm.out, newdata = data.frame(x = x), se.fit = TRUE,interval = 'confidence'))
#   # observeddata[2]
#   dnorm(x =y, mean = predict(lm.out,newdata = data.frame(x = x)),sd = (observeddata[2]))
#   }
fxy2 <- function(x,y){  dnorm(x =y, mean = predict(lm.out,newdata = data.frame(x = x)),sd = (slm$sigma))}
# fxy2 <- function(x,y){  dnorm(x =y, mean = predict(lm.out,newdata = data.frame(x = x)),sd = (sqrt(deviance(lm.out)/df.residual(lm.out))))}
z <- outer(X = xvals,Y = yvals,FUN = fxy2)
z2 <- outer(X = xvals,Y = yvals,FUN = fxy)
# axx <- list(nticks = 4,range = c(-25,75))
# axy <- list(nticks = 4,range = c(-25,75))
# axz <- list(nticks = 4,range = c(0,1))
df.lm$yhat = fitted(lm.out)
# df.lm$yhat = predict(lm.out)
df.lm$z = 0
lm.model.interactive = plotly::plot_ly(x = xrp,
                                       y = yvals,
                                       z = t(z), 
                                       type = "surface",
                                       colors = terrain.colors(100),#c("darkblue", "yellow", "darkred"), 
                                       opacity = 1,
                                       showlegend = F,
                                       width = 500, height = 500) %>% 
  hide_colorbar() %>% 
  add_markers(data = df.lm, x = ~x, y = ~y, z = 0,
              mode = "markers", 
              type = "scatter3d", 
              marker = list(size = 5, color = "black"), 
              name="Points",showlegend = F) %>% 
  add_trace(data = df.lm, x = ~x, y = ~yhat, z = 0, 
            name="Linear", 
            line = list(dash="solid", color = "red", width = 10.5),
            inherit = FALSE,showlegend = F,
            mode = "lines",
            type = "scatter3d") %>% 
  layout(#title = '<b> Linear model expected values </b>',
         # margin = list(t = 80),
    autosize = F,
    margin =list(l = 0, r=0,b = 0, t = 0),
         scene = list(#xaxis=axx,yaxis=axy,
           xaxis = list(#title = 'X',
             range = c(range(xrp)), autorange = "reversed"),
           yaxis = list(#title = 'Y',
             range = range(yvals), autorange = "reversed"),
           zaxis=list(title = 'Z',nticks = 10,range = c(0,1)),
           camera = list(eye = list(x = 1.8, # 0  see https://plotly.com/python/3d-camera-controls/
                                    y = 1.8, # 0.9
                                    z = 0.5))),#2.0))),
           # camera = list(eye = list(x = 0, # 0  see https://plotly.com/python/3d-camera-controls/
           #                          y = .9, # 0.9
           #                          z = 2))),#2.0))),
         showlegend = F
         )

# lm.model.interactive
## save the output 
# htmlwidgets::saveWidget(lm.model.interactive, file = "images/lm.model.interactive.html")
# path_to_python <- "~/miniconda3/bin/python"
# suppressMessages(reticulate::use_python(path_to_python))
# suppressMessages(reticulate::import("plotly"))
# suppressMessages(reticulate::import("kaleido"))
# plotly::save_image(lm.model.interactive, file = "images/lm.model.interactive.png")


## ----simulate.gaussian.model, echo=FALSE, fig.width=8,fig.height=6--------------------------------------------------------------------------------------------

sim.lm = simulate(lm.out, nsim = 1e5, seed = 12) %>% data.frame
lower_ci_sim <- apply(sim.lm, 1, function(x) quantile(x, probs = 0.025) )
upper_ci_sim <- apply(sim.lm, 1, function(x) quantile(x, probs = 0.975) )
sims_summary <- data.frame(
  lower = lower_ci_sim,
  upper = upper_ci_sim
)

# Plot the data 
plot(y~x, data = df.lm, pch = 19, col = b.5)

# Add permutations
permutate.df =replicate(n = 200, # nperm 
          expr = data.frame(y = sample(df.lm$y,size = nrow(df.lm), replace = FALSE), x = df.lm$x),
          simplify = FALSE)
lm.out.perm = mapply(lm, permutate.df)
invisible(apply(lm.out.perm,2,function(x) abline(x,col = scales::alpha("orange",.5))))
lm.out.perm = mapply(lm, permutate.df)

# Add bootstrap
bootstrap.fun <- function(data) {
  tmp.sample = sample(1:nrow(data),size = c(nrow(data)-20),replace = FALSE)
  data.frame(y = data[tmp.sample,"y"], 
             x = data[tmp.sample,"x"])
}
bootstrap.df =replicate(n = 200, # nperm 
                        expr = bootstrap.fun(df.lm),
                        simplify = FALSE)

lm.out.boot = mapply(lm, bootstrap.df)
invisible(apply(lm.out.boot,2,function(x) abline(x,col = scales::alpha("red",.9))))


### Adding the REAL prediction interval 
new.x =  seq(min(df.lm$x),max(df.lm$x), length.out = 100)
p <- predict(lm.out, 
             newdata = data.frame(x = new.x), se.fit=TRUE,     
             interval="prediction")

polygon(x = c(new.x,rev(new.x)),  
        y = c(p$fit[,"lwr"], rev(p$fit[,"upr"])), col = scales::alpha("red",.5),
        density=10)


# Add a line to the plot 
df.lm$sims = sims_summary
o.sim = df.lm[order(df.lm[,1]),]
polygon(x = c(o.sim$x, rev(o.sim$x)),  
        y = c(o.sim$sims$lower,rev(o.sim$sims$upper)), col = scales::alpha("grey50",.5))



# GeneralStandardDev<-sd(lm.out$residuals)
# UpperLine<- lm.out$coefficients[1]+lm.out$coefficients[2]*df.lm$x + GeneralStandardDev
# LowerLine<- lm.out$coefficients[1]+lm.out$coefficients[2]*df.lm$x - GeneralStandardDev
# 
# lines(df.lm$x, UpperLine, col = "blue")
# lines(df.lm$x, LowerLine, col = "blue")
pred = predict(lm.out, interval = "confidence") # 95%  confidence interval
# lines(df.lm$x[order(df.lm$x)], pred[, "lwr"][order(df.lm$x)], col = "blue") 
# lines(df.lm$x[order(df.lm$x)], pred[, "upr"][order(df.lm$x)], col = "blue") 


polygon(x = c(df.lm$x[order(df.lm$x)], rev(df.lm$x[order(df.lm$x)])),  
        y = c(pred[, "lwr"][order(df.lm$x)],rev(pred[, "upr"][order(df.lm$x)])), 
        col = scales::alpha("blue",.5))
abline(lm.out)
points(df.lm$x,df.lm$y)
points(mean(df.lm$x),mean(df.lm$y), col = "red", pch = 19)
abline(v =mean(df.lm$x), h = mean(df.lm$y), lty = 3)



## ----glm_assumtion_expected_mean, echo=FALSE, fig.width=8,fig.height=5----------------------------------------------------------------------------------------
### Linear model 
par(mfrow=c(1,2), mar=c(4,4,2,2), cex = 1.2)
set.seed(12345678)
n = 100; beta0 = 2.5; beta1 = 0.8
x.lm = rnorm(n = n, mean = 10, sd = 1)
err = rnorm(n = n, mean = 0, sd = 1)
# Linear combination 
y.lm = beta0 + beta1*x.lm + err
# Make a dataframe of the data 
df.lm = data.frame(x = x.lm, y = y.lm)
# Colour 
b.5 = scales::alpha("black",alpha = .5)

# PLot the data 
plot(y~x, data = df.lm, pch = 19, col = b.5, xlab = "", cex = 0.5, main = "Linear regression")
title(xlab = "x", ylab="", line=2.2, cex.lab=1.2)
# Model the data 
lm.out = lm(y~x, data = df.lm)
# Add a line to the plot 
abline(lm.out)
####
yvals = seq(floor(range(df.lm$y)[1]),ceiling(range(df.lm$y)[2]),length.out = 1000)
xvals = seq(8,12, length.out = 10)
slm=summary(lm.out)
yh = dnorm(x =yvals, mean = predict(lm.out,newdata = data.frame(x = xvals)),sd = (slm$sigma))
xrp = rep(xvals, length(yh)/length(xvals)) # each=length(yh)/length(xvals)
points(xrp,yvals, cex = 1.5*yh/max(yh), pch = 19, col = scales::alpha("red",.5))

#####
set.seed(42)
n = 500
x = rnorm(n = n, mean = 0, sd = 1)
# Rescale the data
xz = scale(x)
log.mu = 2 + 0.8*xz
y = rpois(n = n, lambda = exp(log.mu)) 

# Combine the data in a dataframe 
df = data.frame(y = y, x = x)

#now feed it to glm:
glm.poisson = glm( y~x, data=df, family="poisson")
plot(y~x, data = df, col = scales::alpha("black", 0.5), pch = 19, cex = 0.5, main = "Poisson regression")
newdata <- data.frame(x = seq(min(x), max(x), len = n))
newdata$y = predict(object = glm.poisson, newdata = newdata, type = "response") 
lines(x = newdata$x,
      y = newdata$y, col = "red",lwd = 2)

yvals = seq(floor(range(df$y)[1]),ceiling(range(df$y)[2]),by = 1)
xvals = seq(-3,3, length.out = 20)
tmp.pois = NULL
pred.lbda= exp(predict(glm.poisson,newdata = data.frame(x = xvals)))
for (i in 1:length(pred.lbda)) {
  yh = dpois(x =yvals, lambda = pred.lbda[i])
  tmp.pois = cbind(tmp.pois,yh)
}
xrp = xvals
for (j in 1:ncol(tmp.pois)) {
  points(rep(xrp[j],length(yvals)),yvals, cex = 1.5*tmp.pois[,j]/max(tmp.pois[,j]), pch = 19, 
         col = scales::alpha(c("red"),.5))# Change the colour depending on the even-odd index
         # col = scales::alpha(c("red", "blue")[j%%2+1],.5))# Change the colour depending on the even-odd index
}



## ----poisson_assumtion_expected_mean_interactive, echo=FALSE, fig.width=3,fig.height=3------------------------------------------------------------------------


library(plotly)
library(dplyr)

beta <- glm.poisson$coefficients
fxy2 <- function(x,y){ dpois(x = y, lambda = exp(predict(glm.poisson, newdata = data.frame(x = xvals))))}
z <- outer(X = xvals,Y = yvals,FUN = fxy2)
df$yhat = fitted(glm.poisson)
df$z = 0
poisson_interactive = plotly::plot_ly(x = xrp,
                y = yvals,
                z = t(z), 
                type = "surface",
                colors = terrain.colors(100),#c("darkblue", "yellow", "darkred"), 
                opacity = 1,
                showlegend = F, width = 500, height = 500) %>% 
  hide_colorbar() %>% 
  add_markers(data = df, x = ~x, y = ~y, z = 0,
              mode = "markers", 
              type = "scatter3d", 
              marker = list(size = 5, color = "black"), 
              name="Points",showlegend = F) %>% 
  add_trace(x = xvals, y = exp(predict(glm.poisson,newdata = data.frame(x = xvals))), z = 0,
            name="Linear",
            line = list(dash="solid", color = "red", width = 10.5),
            inherit = FALSE,showlegend = F,
            mode = "lines",
            type = "scatter3d") %>%
  layout(#title = '<b> Linear model expected values </b>',
    autosize = F,
    margin =list(l = 0, r=0,b = 0, t = 0),
    scene = list(#xaxis=axx,yaxis=axy,
      xaxis = list(#title = 'X',
        range = c(range(xrp)), autorange = "reversed"),
      yaxis = list(#title = 'Y',
        range = range(yvals), autorange = "reversed"),
      zaxis=list(title = 'Z',nticks = 10,range = c(0,1)),
      camera = list(eye = list(x = 1.8, # 0  see https://plotly.com/python/3d-camera-controls/
                               y = 1.8, # 0.9
                               z = 0.7))),#2.0))),
    # camera = list(eye = list(x = 0, # 0  see https://plotly.com/python/3d-camera-controls/
    #                          y = .9, # 0.9
    #                          z = 2))),#2.0))),
    showlegend = F
  )
# poisson_interactive
# # save the output
# htmlwidgets::saveWidget(poisson_interactive, file = "images/poisson_interactive.html")
# path_to_python <- "~/miniconda3/bin/python"
# suppressMessages(reticulate::use_python(path_to_python))
# suppressMessages(reticulate::import("plotly"))
# suppressMessages(reticulate::import("kaleido"))
# plotly::save_image(poisson_interactive, file = "images/poisson_interactive.png")


## ----Sim_t_test_anova-----------------------------------------------------------------------------------------------------------------------------------------
set.seed(1234); n = 1000
y1 = rnorm(n, mean = 15, sd = 1)
y2 = rnorm(n, mean = 15.5, sd = 1)

sim.aov1 = data.frame(y = y1, gr = "A")
sim.aov2 = data.frame(y = y2, gr = "B")
df.aov = rbind(sim.aov1, sim.aov2)
df.aov$gr = factor(df.aov$gr)

# t.test(y~gr, data = df.aov) or
aov.out = aov(y~gr, data = df.aov)
#summary(aov.out)
tk.test = TukeyHSD(aov.out)
round(tk.test$gr,2)


## ----Sim_t_test_anova_plot, fig.width=5,fig.height=5----------------------------------------------------------------------------------------------------------
plot(y~gr, data = df.aov)


## ----simulate_lm_data, echo=-1, fig.width=4,fig.height=2------------------------------------------------------------------------------------------------------
par(mar=c(4,4,0.1,0.1))
set.seed(12345678)
n = 100; beta0 = 2.5; beta1 = 0.8
x.lm = runif(n = n) # or it could be something else! rnorm(n = n, mean = 10, sd = 1)
err = rnorm(n = n, mean = 0, sd = 1)
# Linear combination 
mu = beta0 + beta1*x.lm 
y.lm = mu + err # same as y.lm = rnorm(n = n, mu, sd = 1)
# Make a dataframe of the data 
df.lm = data.frame(x = x.lm, y = y.lm)


## ----simulate_lm_plot, echo=FALSE, fig.width=4,fig.height=3---------------------------------------------------------------------------------------------------
par(mar = c(4,4,.5,.5))
# Colour 
b.5 = scales::alpha("black",alpha = .5)

# PLot the data 
plot(y~x, data = df.lm, pch = 19, col = b.5)

# Model the data 
lm.out = lm(y~x, data = df.lm)
# Add a line to the plot 
abline(lm.out)


## ----simulate_lm_output---------------------------------------------------------------------------------------------------------------------------------------
round(coefficients(summary(lm.out)), 4)
# Adjusted R^2
summary(lm.out)$adj.r.squared
# summary(lm.out)$fstatistic
# anova(lm.out)


## ----simulate_lm_resimulate_from_model, echo=-1, fig.width=13,fig.height=3------------------------------------------------------------------------------------
par(mfrow = c(1,2), mar = c(3.5,4,1,1), cex = 1.2)
sim.lm = simulate(lm.out, nsim = 2000, seed = 12) # simulate based on predicted value and sigma of model #<<
rx = range(c(rowMeans(sim.lm), fitted(lm.out)))
hist(rowMeans(sim.lm),xlim=rx,main="Hist simulation",xlab="");title(xlab="rowMeans(sim.lm)",line=2.2,cex.lab=1.2)
hist(fitted(lm.out), xlim = rx,main="Hist fitted",xlab="");title(xlab="fitted(lm.out)",line=2.2,cex.lab=1.2)
c(mean(rowMeans(sim.lm)), mean(y.lm),  mean(fitted(lm.out))) # compare
rbind(head(rowMeans(sim.lm)), head(fitted(lm.out))) # Average of all simulations (in row) gives the fitted values


## ----lm_simulate_slope, echo=FALSE, fig.width=12,fig.height=4-------------------------------------------------------------------------------------------------
par(mar=c(4,4,1,1), cex = 1.3)

set.seed(12345678)
lm.sim.fun <- function(n = 100, 
                       mean = 10, sd = 1, sd.err = 1,
                       beta0 = 2.5, beta1 = 0.8,
                       plot = FALSE, ret = TRUE) {
  # Generate x values 
  x.lm = rnorm(n = n, mean = mean, sd = sd)
  # Generate error for each x value 
  err = rnorm(n = n, mean = 0, sd = sd.err)
  
  # Linear combination 
  y.lm = beta0 + beta1*x.lm + err
  
  # Make a dataframe of the data 
  df.lm = data.frame(x = x.lm, y = y.lm)
  
  # Colour 
  b.5 = scales::alpha("black",alpha = .5)
  
  # Model the data 
  lm.out = lm(y~x, data = df.lm)
  # Get the summary 
  s.out = summary(lm.out)
  
  # Might want to add the plot 
  if (plot) {
    # PLot the data 
    plot(y~x, data = df.lm, pch = 19, col = b.5)
    # Add a line to the plot 
    abline(lm.out)
  } # end if plot
  
  # Return some interesting parameters for future analysis
  if (ret) {
  return(list(lm.out = s.out,
              int = s.out$coefficients["(Intercept)",1],
              slope = s.out$coefficients["x",1]))
  }
}

# Number of repetitions to replicate our simulation
nb.rep = 1000
# Make an empty data object to record what we are interested in 
df.all.sim = NULL
beta1 = 0.8
for (i in c(5,20,200)) {
  l.rp = replicate(n = nb.rep, 
                   expr = lm.sim.fun(n = i,beta1 = beta1), 
                   simplify = FALSE)
  all.int = unlist(lapply(l.rp, function(x) x$int)) # Get intercept
  all.slp = unlist(lapply(l.rp, function(x) x$slope)) # Get slope
  df.all.sim = rbind(df.all.sim,
                     data.frame(n = i, # number of replicate 
                                all.int, # all intercepts
                                all.slp)) # all slopes 
}

plot.dens <- function(density, data, subset, main) {
  plot(density,   xlim = c(-3,4), main = main, cex = 1.3)
abline(v = mean(data[data$n==subset,"all.slp"]),   lty =1)
abline(v = beta1, lty = 3)
polygon(x = c(-10, density$x[density$x>-10 & density$x < 10], 10), 
        y = c(0, density$y[density$x>=-10 & density$x <= 10], 0),
        col=scales::alpha("blue",.5))
}

par(mfrow = c(1,3), cex = 1.2)
# get the density 
dens5 = density(df.all.sim[df.all.sim$n==5,"all.slp"]  )
dens20 = density(df.all.sim[df.all.sim$n==20,"all.slp"] )
dens200 = density(df.all.sim[df.all.sim$n==200,"all.slp"])

# Plot for each parameter we are interested in 
plot.dens(density = dens5, data = df.all.sim,subset = 5, main = paste("slope =",beta1,"n = 5"))
plot.dens(density = dens20, data = df.all.sim,subset = 20, main = paste("slope =",beta1,"n = 20"))
plot.dens(density = dens200, data = df.all.sim,subset = 200, main = paste("slope =",beta1,"n = 200"))

# lm.f <- function(n,mean,sd,b0,b1, plot = FALSE) {
#   x = rnorm(n, mean = mean, sd = sd)
#   beta0 = b0; beta1 = b1
#   y = beta0 + beta1 * x + rnorm(n)
#   lm.out = lm(y~x)
#   s.lm.out = summary(lm.out)
#   if (plot) {
#     plot(x,y, pch = 19)
#     abline(lm.out)
#   }
#   return(s.lm.out)
# }
# 
# par(mfrow = c(1,3))
# for (i in c(5,20,1000)) {
#   rp.lm = replicate(n = 1000, 
#                     expr = lm.f(n = i,mean = 19, sd = 1, b0 = 1, b1 = 4, plot = FALSE)$coefficients[c("(Intercept)","x"),"Estimate"],
#                     simplify = FALSE)
#   sim.est.lm = do.call(rbind, rp.lm)
#   dens = density(sim.est.lm[,"x"])
#   plot(dens, main = paste("n=",i), xlim = c(-1,10))
#   polygon(c(-10, dens$x[dens$x>-10 & dens$x < 10], 10), 
#           c(0, dens$y[dens$x>=-10 & dens$x <= 10], 0), 
#           col=scales::alpha("blue",.5))
# }


## ----lm_simulate_r_sq, echo=FALSE, fig.width=10,fig.height=5--------------------------------------------------------------------------------------------------
par(mfrow = c(2,3),mar=c(4,4,1,1), cex = 1.1)
set.seed(2345)
nb.rep = 1000
sd.err = 2
df.all.sim = NULL
for (i in c(5,20,200)) {
  l.rp = replicate(nb.rep, lm.sim.fun(n = i,beta1 = 2,sd.err = sd.err), simplify = FALSE)
  all.adjr2 = unlist(lapply(l.rp, function(x) x$lm.out$adj.r.squared))
  all.int = unlist(lapply(l.rp, function(x) x$int))
  all.slp = unlist(lapply(l.rp, function(x) x$slope))
  df.all.sim = rbind(df.all.sim,data.frame(n = i, all.int, all.slp,all.adjr2))
}
dens5 = density(df.all.sim[df.all.sim$n==5,"all.adjr2"]  )
dens20 = density(df.all.sim[df.all.sim$n==20,"all.adjr2"] )
dens200 = density(df.all.sim[df.all.sim$n==200,"all.adjr2"])
xlim = c(-1,1.3)
plot(dens5,   xlim = xlim, main = "n = 5");   abline(v = mean(df.all.sim[df.all.sim$n==5,"all.adjr2"]),   lty =1)
polygon(c(-10, dens5$x[dens5$x>-10 & dens5$x < 10], 10), c(0, dens5$y[dens5$x>=-10 & dens5$x <= 10], 0), col=scales::alpha("blue",.5))
plot(dens20,  xlim = xlim, main = "n = 20");  abline(v = mean(df.all.sim[df.all.sim$n==20,"all.adjr2"]),  lty =1)
polygon(c(-10, dens20$x[dens20$x>-10 & dens20$x < 10], 10), c(0, dens20$y[dens20$x>=-10 & dens20$x <= 10], 0), col=scales::alpha("blue",.5))
plot(dens200, xlim = xlim, main = "n = 200"); abline(v = mean(df.all.sim[df.all.sim$n==200,"all.adjr2"]), lty =1)
polygon(c(-10, dens200$x[dens200$x>-10 & dens200$x < 10], 10), c(0, dens200$y[dens200$x>=-10 & dens200$x <= 10], 0), col=scales::alpha("blue",.5))
lm.sim.fun(n = 5,  sd.err = sd.err,plot = T, ret = F)
lm.sim.fun(n = 20, sd.err = sd.err,plot = T, ret = F)
lm.sim.fun(n = 200,sd.err = sd.err,plot = T, ret = F)



## ----lm_simulate_2x, echo=FALSE, fig.width=8,fig.height=4-----------------------------------------------------------------------------------------------------
set.seed(12345678)
n = 100
x1 = rnorm(n = n, mean = 10, sd = 1)
x2 = rnorm(n = n, mean = 15, sd = 6)
# x1 = scale(x1)
# x2 = scale(x2)
beta0 = 2.5
beta1 = .8
beta2 = .5
err = rnorm(n = n, mean = 0, sd = 1)

# Linear combination 
y.lm = beta0 + beta1*x1 + beta2*x2 + err

# Make a dataframe of the data 
df.lm = data.frame(x1 = x1, x2 = x2, y = y.lm)

# Colour 
b.5 = scales::alpha("black",alpha = .5)

par(mfrow=c(1,2))
# PLot the data 
plot(y~x1, data = df.lm, pch = 19, col = b.5)

# Model the data 
lm.out = lm(y~x1+x2, data = df.lm)
# Add a line to the plot 
abline(a = coef(lm.out)["(Intercept)"]+mean(x2)*coef(lm.out)["x2"], b = coef(lm.out)["x1"])

# PLot the data 
plot(y~x2, data = df.lm, pch = 19, col = b.5)
# Add a line to the plot 
abline(a = coef(lm.out)["(Intercept)"]+mean(x1)*coef(lm.out)["x1"], b = coef(lm.out)["x2"])
# abline(lm.out)


## ----lm_simulate_2x_out, echo=FALSE---------------------------------------------------------------------------------------------------------------------------
cat("Coefficients table:\n"); round(coefficients(summary(lm.out)), 4)
# Adjusted R^2
cat("R Squared: "); summary(lm.out)$adj.r.squared
# summary(lm.out)$fstatistic
# anova(lm.out)
# vcov(lm.out)


## ----sim_lm_cat, fig.width=5,fig.height=5---------------------------------------------------------------------------------------------------------------------
set.seed(12345678)
n = 100; meanx = 10
x.lm = rnorm(n = n, mean = meanx, sd = 1)
cat.lm = c("Control", "Treatment")
gr = gl(n = length(cat.lm), k = n/2, labels = cat.lm) # makes factors (default)
# Make a dummy matrix
dummy.mat = matrix(data = 0,
                   nrow = n,
                   ncol = length(cat.lm))
colnames(dummy.mat) <- paste0(c(cat.lm))
dummy.mat[,1] <- 1 # Intercept
dummy.mat[,2] <- ifelse(gr=="Treatment",1,0)

beta1 = 0.8; betaControl = .2; betaTreament = .7
# Combine the data to get the linear prediction 
lin.pred.x = dummy.mat %*% c(betaControl,betaTreament) + beta1 * x.lm
y.lm = lin.pred.x + rnorm(n = n, mean = 0, sd = 1)
df.lm = data.frame(x = x.lm, y = y.lm, gr = gr)
b.5 = scales::alpha("black",alpha = .5); r.5 = scales::alpha("red",alpha = .5)
lm.out = lm(y~x+gr, data = df.lm);   lmcoef <- lm.out$coefficients


## ----sim_lm_cat_out, echo=FALSE, eval=TRUE, results = 'hide'--------------------------------------------------------------------------------------------------
# Think of beta1 as multiplying the mean of X by beta1 and then seeing the effect of the betaControl and betaTreament
# See this 
mean(x.lm) # Mean of the X variable 
mean(y.lm) # Mean of the Y variable 

# Theoretical
t.x = meanx * beta1
# Observed 
o.x = mean(x.lm) * lm.out$coefficients["x"]

# Theoretical
t.xc = meanx * beta1 + betaControl
# Observed 
o.xc = mean(x.lm) * lm.out$coefficients["x"] + lm.out$coefficients["(Intercept)"]

# Theoretical
t.xt = meanx * beta1 + betaControl + betaTreament
# Observed 
o.xt = mean(x.lm) * lm.out$coefficients["x"] + lm.out$coefficients["(Intercept)"] + lm.out$coefficients["grTreatment"]
mat.exp.obs = matrix(c(meanx, mean(x.lm),
                       NA,NA,
                       t.x,o.x,
                       t.xc,o.xc,
                       t.xt,o.xt,
                       beta1,lm.out$coefficients["x"],
                       betaControl, lm.out$coefficients["(Intercept)"],
                       betaTreament, lm.out$coefficients["grTreatment"]
                       ), ncol = 2, byrow = TRUE, 
                     dimnames = list(c("x","y",
                                       "m.beta1","m.beta1.bC","m.beta1.bC.bT", 
                                       "beta1","betaC","betaT"),
                                     c("Expected","Observed")))
# m.beta1 # is the global intercept
df.eo = as.data.frame(mat.exp.obs)
df.eo$Expected[which(is.na(df.eo$Expected))] <- list(unique(dummy.mat %*% c(betaControl,betaTreament) + beta1 * meanx))
df.eo$Observed[which(is.na(df.eo$Observed))] <- list(aggregate(x = df.lm[,c("y")],by = list(df.lm[,c("gr")]), mean)$x)

df.eo$Notes = c("Mean of X","Mean of each category in Y", "Global intercept (mean)","Global + control","Global + control + treatment","","","")
# Here you can see the beta1 (x), betaControl (Intercept) and betaTreament (grTreatment)
summary(lm(y~x+gr, data = df.lm))

# Get JUST the DIFFERENCE with the treatment (Intercept is the Control)
gr.means = lm(y~gr, data = df.lm) 
summary(gr.means)

# Calculate the mean of the 2 groups 
summary(lm(y~gr-1, data = df.lm) )


# Another way to calculate the mean of the 2 groups 
aggregate(x = df.lm[,c("y")],by = list(df.lm[,c("gr")]), mean)

# Compare to this (ANCOVA!)
summary(aov(y~x+gr,data = df.lm))
anova(lm.out)


## ----sim_lm_cat_out_summary-----------------------------------------------------------------------------------------------------------------------------------
# Here you can see the beta1 (x), betaControl (Intercept) and betaTreament (grTreatment)
summary(lm(y~x+gr, data = df.lm))
# Get JUST the DIFFERENCE with the treatment (Intercept is the Control)
# gr.means = lm(y~gr, data = df.lm) 
# coefficients(gr.means)


## ----sim_lm_cat_plot, echo=-1, fig.width=8,fig.height=3.5-----------------------------------------------------------------------------------------------------
par(mfrow=c(1,2), mar = c(4,4,.5,.5))
plot(y~x, data = df.lm, pch = 19, col = ifelse(gr=="Treatment",r.5,b.5))
abline(a= lmcoef["(Intercept)"], b= lmcoef["x"], col = "black") # control
abline(a= lmcoef["(Intercept)"] + lmcoef["grTreatment"], b=lmcoef["x"], col="red")
boxplot(y~gr, data = df.lm, pch = 19, col = ifelse(gr=="Treatment",r.5,b.5))
round(coefficients(summary(lm.out)), 4) # betaControl is (Intercept)
summary(lm.out)$adj.r.squared           # betaTreament is grTreatment. 


## ----show_df, echo=FALSE--------------------------------------------------------------------------------------------------------------------------------------
df.eo


## ----sim_lm_interaction_discrete, echo=FALSE, fig.width=8,fig.height=3----------------------------------------------------------------------------------------
set.seed(1)
n = 100
x.lm1 = rbinom(n = n, size = 1, prob = 0.5)
x.lm2 = rnorm(n = n, mean = 1, sd = 1) # Note that this doesn't have to be normally distributed. This could be a uniform distribution or from a binomial.

# x.lm2 = scale(x.lm2)
beta0 = 2.5
beta1 = 1.5
beta2 = 2
beta3 = -3
err.lm = rnorm(n = n, mean = 0, sd = 1)
y.lm = beta0 + beta1*x.lm1 + beta2*x.lm2 + beta3*x.lm1*x.lm2 + err.lm


df.lm = data.frame(x1 = x.lm1, x2 = x.lm2, y = y.lm)
b.5 = scales::alpha("black",alpha = .5)
par(mfrow=c(1,3))
# apply(df.lm, 2, hist)
hist(df.lm$x1)
hist(df.lm$x2)
hist(df.lm$y)

lm.out = lm(y~x1*x2, data = df.lm)
# summary(lm.out)
round(coefficients(summary(lm.out)), 4)
# Adjusted R^2
summary(lm.out)$adj.r.squared
# summary(lm.out)$fstatistic
# anova(lm.out)
# coef(lm.out)


## ----lm_sim_interaction, echo=-1, fig.width=5,fig.height=5----------------------------------------------------------------------------------------------------
par(mfrow=c(1,1), mar = c(4,4,.1,.1))
plot(x = df.lm[df.lm$x1 == 0, ]$x2, y = df.lm[df.lm$x1 == 0, ]$y, # Add x1 = 0
     col = rgb(red = 0, green = 0, blue = 1, alpha = 0.25), pch = 19,
     xlab = "x2", ylab = "y", ylim = range(df.lm$y)); abline(h=0, v=0, lt =3)
abline(a = coef(lm.out)[1], b = coef(lm.out)[3], col = "blue", pch = 19, lwd =2)
points(x = df.lm[df.lm$x1 == 1, ]$x2, y = df.lm[df.lm$x1 == 1, ]$y, # Add x1 = 1
       col = rgb(red = 1, green = 0, blue = 0, alpha = 0.25), pch = 19)
abline(a = coef(lm.out)[1] + coef(lm.out)[2], 
       b = coef(lm.out)[3] + coef(lm.out)[4], col = "red", lwd = 2)


## ----sim_lm_interaction, echo=FALSE, fig.width=8,fig.height=5-------------------------------------------------------------------------------------------------
library(viridis)
set.seed(1)
n = 100
# x.lm1 = rbinom(n = n, size = 1, prob = 0.5)
x.lm1 = rnorm(n = n, mean = 5, sd = 1)
x.lm2 = rnorm(n = n, mean = 2, sd = 1) # Note that this doesn't have to be normally distributed. This could be a uniform distribution or from a binomial.
# x.lm1 = scale(x.lm1)
# x.lm2 = scale(x.lm2)
beta0 = 2.5
beta1 = 1.5
beta2 = 2
beta3 = 3
err.lm = rnorm(n = n, mean = 0, sd = 1)
y.lm = beta0 + beta1*x.lm1 + beta2*x.lm2 + beta3*x.lm1*x.lm2 + err.lm

df.lm = data.frame(x1 = x.lm1, x2 = x.lm2, y = y.lm)
b.5 = scales::alpha("black",alpha = .5)


## ---- echo=-c(1:2), fig.width=8,fig.height=2,.5---------------------------------------------------------------------------------------------------------------
par(mfrow=c(1,3))
# apply(df.lm, 2, hist)
hist(df.lm$x1);hist(df.lm$x2);hist(df.lm$y)
lm.out = lm(y~x1*x2, data = df.lm); lmcoef = coef(lm.out)
# Make a new range of x2 values on which we will test the effect of x1 
x2r = range(x.lm2); x2.sim = seq(x2r[1], x2r[2], by = .5); x2s.lg=length(x2.sim)
# This is the effect of x1 at different values of x2 (moderates the effect of x1)
eff.x1 <- lmcoef["x1"] + lmcoef["x1:x2"] * x2.sim # gets slopes  
eff.x1.int <- lmcoef["(Intercept)"] + lmcoef["x2"] * x2.sim # gets  int.  
eff.dat <- data.frame(x2.sim, eff.x1, eff.x1.int) # Put that in dataframe 
virPal <- viridis::viridis(length(x2.sim), alpha = .8) # Get colours 
eff.dat$col <- virPal[as.numeric(cut(eff.dat$x2.sim, breaks = x2s.lg))]
df.lm$col <- virPal[as.numeric(cut(df.lm$x2, breaks = x2s.lg))]
df.lm$line <- c("black","red")[as.numeric(cut(df.lm$x2, breaks = x2s.lg)) %% 2+1]


## ----lm_sim_interaction_plot_1x, echo=-c(1:3), fig.width=8,fig.height=4---------------------------------------------------------------------------------------
# plot(x = eff.dat$x2.sim, y = eff.dat$eff.x1, type = "l",#col =df.lm$line,
#      # xlim=range(x2.sim),
#      pch = 19, xlab = "Level of x2", ylab = "Marginal effect of x1")
par(mfrow=c(1,1), mar =c(4,4,1,1))
plot(x = df.lm$x1, y = df.lm$y, bg = df.lm$col, 
     pch = 21, xlab = "x1", ylab = "y")
apply(eff.dat, 1, function(x) abline(a = x[3], b = x[2], col = x[4], lwd  = 2))
abline(h = 0, v = 0,lty = 3)
legend("topleft", title = "x2",legend = round(eff.dat$x2.sim,1), lty = 1, lwd = 3,
       col = eff.dat$col, bg = scales::alpha("white",.5))


## ----lm_sim_interaction_out_summary---------------------------------------------------------------------------------------------------------------------------
summary(lm.out)


## ----challenge4.one.rep, fig.width=5,fig.height=5-------------------------------------------------------------------------------------------------------------
set.seed(1234)
n = 10;  sd.e = 0.3
control = rnorm(n, mean = 5.03, sd = 0.28)
treatmt = rnorm(n, mean = 4.66, sd = 0.49)
gr <- gl(n = 2, k = n, length = n*2, labels = c("ctl","trt"))
weight = c(control,treatmt) + rnorm(2*n, mean = 0, sd = sd.e)
plant.weight = data.frame(weight=weight, gr = gr)
lm.out.int <- lm(weight ~ gr, data = plant.weight)
lm.out.noi <- lm(weight ~ gr-1, data = plant.weight) # Comparing the data to 0 (int. = 0, which becomes the reference)
anova(lm.out.int) # this time, it was significant 
# summary(lm.out.noi) # Just to see the estimate (NOT THE P-VALUES)
s.lm.out = summary(lm.out.int)
# You can get the p-value with this function 
pt(q = abs(s.lm.out$coefficients["grtrt","t value"]),df = s.lm.out$df[2], lower.tail = F) *2


## ----challenge4.diagnose1, echo=-c(1), fig.width=5,fig.height=4-----------------------------------------------------------------------------------------------
par(mfrow = c(1,1), mar = c(4, 4, .1, 0.1))
plot(weight~gr, data=plant.weight, las=1)


## ----challenge4.diagnose2, echo=-c(1), fig.width=5,fig.height=5-----------------------------------------------------------------------------------------------
par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
plot(lm.out.int, las = 1)


## ----challenge4.fun, fig.width=8,fig.height=5-----------------------------------------------------------------------------------------------------------------
# Make the function for the simulation
exp.plant <- function(n = 100, # number of seeds in EACH group (2*n = total)
                      sd.e = 0.3, plot=F, ret = T) {
  control = rnorm(n, mean = 5.03, sd = 0.28)
  treatmt = rnorm(n, mean = 4.66, sd = 0.49)
  gr <- gl(n = 2, k = n, length = 2*n, labels = c("ctl","trt"))
  weight = c(control,treatmt) + rnorm(2*n, mean = 0, sd = sd.e)
  plant.weight = data.frame(weight=weight, gr = gr)
  if (plot) {
    plot(weight ~ gr, data = plant.weight, las = 1)
  }
  lm.out.int <- lm(weight ~ gr, data = plant.weight)
  s.lm.out = summary(lm.out.int)
  if (ret) {
   return(s.lm.out) 
  }
}



## ----challenge4.fun.plot, echo=-c(1:3), fig.width=8,fig.height=3.5--------------------------------------------------------------------------------------------
par(mfrow = c(2,2), mar = c(4,4,.5,.5))
layout.matrix <- matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2)
layout(layout.matrix, widths = c(2,2), heights = c(1,1))
set.seed(2345); nb.rep = 2000 # number of replications in the simulation 
l.rp = replicate(n = nb.rep, simplify = FALSE,
                 expr = exp.plant( n = 10)$coefficients["grtrt","t value"])
p.val.lm = pt(q = abs(unlist(l.rp)),df = s.lm.out$df[2], lower.tail = F)*2 # Get p-values 
exp.plant(n = 10, plot = T, ret = F) # plot a simulation example 
hist(unlist(l.rp), main = "t-values", probability = T, xlab ="")
lines(density(unlist(l.rp)), col="blue", lwd=2); abline(v = qt(0.025, df = s.lm.out$df[2]))
hist(p.val.lm, main = "p-values", probability = T, xlab ="", xlim = c(0,1))
lines(density(p.val.lm), col="red", lwd=2); abline(v = 0.05)


## ----chall.6.output-------------------------------------------------------------------------------------------------------------------------------------------
qt(0.05, df = s.lm.out$df[2])
sum(unlist(l.rp)<qt(0.025, df = s.lm.out$df[2]))/length(unlist(l.rp))
sum(p.val.lm<=0.05)/length(p.val.lm)


## ----sim_logistic---------------------------------------------------------------------------------------------------------------------------------------------
set.seed(987654)
n = 100
x1 = rnorm(n = n, mean = 6, sd = 1)
x2 = rnorm(n = n, mean = 0, sd = 1)
# Rescale the data
x1z = scale(x1)
x2z = scale(x2)
z = 0 + 2*x1z + 3*x2 # This gives the LOG odds. Recall that $p = odds/(1+odds)$
pr = 1/(1+exp(-z)) # Transform to get the LOG(odds) # inverse-logit function; 
y = rbinom(n = n, size = 1, prob = pr) # Bernoulli response variable 

# Combine the data in a dataframe 
df = data.frame(y = y, x1 = x1, x2 = x2)

# Compute the model
glm.logist = glm( y~x1+x2, data=df, family="binomial")



## -------------------------------------------------------------------------------------------------------------------------------------------------------------
# source : https://sebastiansauer.github.io/convert_logit2prob/
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}
coef(glm.logist)
round(logit2prob(coef(glm.logist)),5)


## ----sim_logistic_glm, echo=-1, fig.width=5,fig.height=3------------------------------------------------------------------------------------------------------
par(mar = c(4,4,0.2,0.1))
plot(y~x1, data = df, col = scales::alpha("black",.5), pch = 19)
newdata <- data.frame(x1=seq(min(x1), max(x1),len=n), 
                      x2 = seq(min(x2), max(x2),len=n))
newdata$y = predict(object = glm.logist, newdata = newdata, type = "response") 
lines(x = newdata$x1, y = newdata$y, col = "red",lwd = 2)


## ----glm_logistic_summary, echo=FALSE-------------------------------------------------------------------------------------------------------------------------
glm.sum = summary(glm.logist) # Note that the output of the GLM is in "logit"
coefficients(glm.sum)[,"Estimate"]
# glm.sum$coefficients[,1]/glm.sum$coefficients[,2] # Estimate/Std. Error = z val.


## ----glm_logistic_odd_prob_linear, echo=FALSE, fig.width=8,fig.height=3---------------------------------------------------------------------------------------
set.seed(6)
n = 200
x1 = rnorm(n = n, mean = 6, sd = 3)
# Rescale the data
x1z = scale(x1)
beta1 = 2
z = 0 + beta1*x1z  
pr = 1/(1+exp(-z))
y = rbinom(n = n, size = 1, prob = pr) 

# Combine the data in a dataframe 
df = data.frame(y = y, x1 = x1)
model.logistic = data.frame(x1,pr,y)
model.logistic = model.logistic[order(model.logistic$x1),]

#now feed it to glm:
glm.logist = glm( y~x1, data=df, family="binomial")
glm.sum = summary(glm.logist)
# cat("Coefficients\n")
glm.sum$coefficients

par(mfrow=c(1,3))
b.5 = scales::alpha("black",.5)
plot(z~x1z, ylab = "Log Odds", main = "Theoretical logistic in log(odds)",
     pch = 19, col = b.5, xlim = c(-5,5), ylim = c(-12,11), asp = 1)
abline(a = 0,
       b = beta1, col = "red")
abline(h=0, v=0,lty = 3)

newdata <- data.frame(x1=seq(min(x1), max(x1),len=n))
newdata$y = predict(object = glm.logist, newdata = newdata, type = "response") 
# log(newdata$y/(1-newdata$y))
plot(z~x1, ylab = "Log Odds", main = "Estimated logistic in log(odds)",
     pch = 19, 
     col = scales::alpha(c("blue","red")[y+1],.3),
     # col = b.5,
     xlim = c(0,10), ylim = c(-12,11), asp = 1)
abline(h=0, v=0,lty = 3)
# newdata$y/(1-newdata$y)
points(x = newdata$x1, 
       y=log(newdata$y/(1-newdata$y)), 
       pch = 19, 
       col = scales::alpha("green",.2))
       # col = c("blue","red")[model.logistic$y+1])
abline(a = glm.sum$coefficients[1,1],
       b = glm.sum$coefficients[2,1])
points(x = 0, y=glm.sum$coefficients[1,1], pch = 19, col = "red")
text(x = 0, y=glm.sum$coefficients[1,1], labels = c("Intercept"), pos =4)

# par(mfrow = c(1,1))
plot(y~x1, data = df, main = "Estimated logistic in probability",
     ylab = "Probability of outcome",
     col = scales::alpha(c("blue","red")[y+1],.5), pch = 19)
abline(h=0.5, v=mean(x1),lty = 3)
lines(x = newdata$x1,
      y = newdata$y, col = "black",lwd = 2)

lines(x = model.logistic$x1,
      y = model.logistic$pr, col = "red",lwd = 2)

myx1 = seq(min(x1),max(x1), length.out = length(x1))
# Rescale the data
my.x1z = scale(myx1)
beta1 = 2
my.z = 0 + beta1*my.x1z
new.pr = 1/(1+exp(-my.z))

new.p = dbinom(c(1,0),size = 1,prob = new.pr)
p = dbinom(c(1,0),size = 1,prob = pr)
new.q = 1-new.p
q = 1-p
df = data.frame(y = y, x1 = x1, x1z = x1z, p, q, myx1, new.pr,pr, new.p)
# points(x = df$x1,
#        y = df$p, cex = 2*df$p, col = scales::alpha(c("green","orange"),.5), pch = 19)

my.x=seq(range(x1z)[1],range(x1z)[2], length.out = 20)

my.pr = 1/(1+exp(-(0 + beta1*my.x)))
p = dbinom(c(1,0),size = 1,prob = rep(my.pr,each =2 ))
q = 1-p
df = data.frame( p, q,my.x= rep(my.x, each =2),my.pr)
points(x = (df$my.x*attr(x1z,which = "scaled:scale"))+attr(x1z,which = "scaled:center"),
       y = df$p, cex = 1.3*exp(df$p), col = scales::alpha(c("green","orange"),.5), pch = 19)
x.vals = (my.x*attr(x1z,which = "scaled:scale"))+attr(x1z,which = "scaled:center")
abline(v = x.vals[seq(1,length(x.vals),by = 2)], lty = 3, col = "grey80")

cat("Probability at intercept is",1/(1+exp(-glm.sum$coefficients[1,1])),"\n")
cat("Probability for an increase in 1 sd of X is",1/(1+exp(-glm.sum$coefficients[2,1])))


## ----glm_logistic_odd_prob_quad, echo=FALSE, fig.width=8,fig.height=3-----------------------------------------------------------------------------------------
set.seed(6)
n = 150
x1 = rnorm(n = n, mean = 6, sd = 3)
# Rescale the data
x1z = scale(x1)
B0 = 1
B1 = .4
B2 = -.5
z = B0 + B1*x1z + B2*x1z^2 # This gives the LOG odds. Recall that $p = odds/(1+odds)$
pr = 1/(1+exp(-z)) # Transform to get the LOG(odds) # inverse-logit function; 
y = rbinom(n = n, size = 1, prob = pr) # Bernoulli response variable 

# Combine the data in a dataframe 
df = data.frame(y = y, x1 = x1z, x = x1)

# Compute the model
glm.logist = glm( y~x1+I(x1^2), data=df, family="binomial")

#now feed it to glm:
glm.sum = summary(glm.logist)
# cat("Coefficients\n")
glm.sum$coefficients

par(mfrow=c(1,3))
b.5 = scales::alpha("black",.5)
r.5 = scales::alpha("red",.5)
col.line = "green"
plot(z~x1z, ylab = "Log Odds", main = "Theoretical logistic in log(odds)",
     pch = 19, col = ifelse(df$y==1,b.5,r.5), xlim = c(-5,5), ylim = c(-12,11)); abline(h=0, v=0,lty = 3)
x.sq = seq(-12,12,length.out = n)
Y = B0 + B1 *x.sq + B2 *x.sq^2
lines(x.sq, Y, col = col.line)
# Add the point of inflection
# points(-B1/(2*B2),c(B0+B1*(-B1/(2*B2))+B2*(-B1/(2*B2))^2), col = "pink", cex = 3, pch = 19)
# Where the log odds =0 BUT without the B0
# points(-B1/(B2),c(B0+B1*(-B1/(B2))+B2*(-B1/(B2))^2), col = "pink", cex = 3, pch = 19)
nx = -B1/B2
B0 + B1 *nx + B2 *nx^2

plot(z~x1, ylab = "Log Odds", main = "Estimated logistic in log(odds)",
     pch = 19, col = ifelse(df$y==1,b.5,r.5), xlim = c(0,10), ylim = c(-12,11));abline(h=0, v=0,lty = 3)

glm.logist2 = glm( y~x+I(x^2), data=df, family="binomial")
glm.sum = summary(glm.logist2)
b0 = glm.sum$coefficients["(Intercept)","Estimate"]
b1 = glm.sum$coefficients["x","Estimate"]
b2 = glm.sum$coefficients["I(x^2)","Estimate"]
x.sq = seq(-2,12,length.out = n)
yhat = b0 + b1 *x.sq + b2 *x.sq^2
lines(x.sq, yhat, col = col.line)
points(x = 0, y=b0, pch = 19, col = "blue")
text(x = 0, y=b0, labels = c("Intercept"), pos =4)

plot(y~x, data = df, main = "Estimated logistic in probability",
     ylab = "Probability of outcome", xlim =range(-5,x1,20),
     col = ifelse(df$y==1,b.5,r.5), pch = 19)
abline(h=0.5, v=mean(x1),lty = 3)
newdata <- data.frame(x=seq(min(-10,x1), max(x1,20),len=n))
newdata$y = predict(object = glm.logist2, newdata = newdata, type = "response") 
lines(x = newdata$x,
      y = newdata$y, col = col.line,lwd = 2)

cat("Probability at intercept is",1/(1+exp(-glm.sum$coefficients[1,1])),"\n")

mod.test <- function(b0,b1,b2,x) {
  b0 + b1 *x + b2 *x^2  
}

# Show the points at certain x values 
points(mean(df$x)+0*sd(df$x),1-1/(1+exp(mod.test(b0,b1,b2,mean(df$x)+0*sd(df$x)))))
points(mean(df$x)+1*sd(df$x),1-1/(1+exp(mod.test(b0,b1,b2,mean(df$x)+1*sd(df$x)))))
points(mean(df$x)+2*sd(df$x),1-1/(1+exp(mod.test(b0,b1,b2,mean(df$x)+2*sd(df$x)))))


## ----glm_logistic_odd_prob_quad_showval-----------------------------------------------------------------------------------------------------------------------
mean(df$x)+0*sd(df$x)
1-1/(1+exp(mod.test(b0,b1,b2,mean(df$x)+0*sd(df$x))))
mean(df$x)+1*sd(df$x)
1-1/(1+exp(mod.test(b0,b1,b2,mean(df$x)+1*sd(df$x))))
mean(df$x)+2*sd(df$x)
1-1/(1+exp(mod.test(b0,b1,b2,mean(df$x)+2*sd(df$x))))
mean(df$x)+5*sd(df$x)
1-1/(1+exp(mod.test(b0,b1,b2,mean(df$x)+5*sd(df$x))))


## ----glm_logistic_cat_1, echo=FALSE---------------------------------------------------------------------------------------------------------------------------
set.seed(12345678)
n <- 10000
# generate categories 
categories = c("A","B", "C")

x <- sample(x = categories, 
            size = n, 
            replace = TRUE, 
            prob = rep(x = 1/length(categories), length(categories)))

# Here is another way to generate factors 
# cat.fact = gl(n = 3,k = 1,length = n,labels = categories)

# Make a dummy matrix
dummy.mat = matrix(data = 0,
                   nrow = length(x),
                   ncol = length(categories))
colnames(dummy.mat) <- paste0("x",c(categories))
dummy.mat[,1] <- 1
dummy.mat[,2] <- ifelse(x=="B",1,0)
dummy.mat[,3] <- ifelse(x=="C",1,0)

# contrasts(as.factor(x))

# Easier
library(dummies)
# dummy(x)

# Set some coefficients
beta0 <- 2
betaB <- 3.1
betaC <- -4.35
beta1 <- 2

xxx = rnorm(n = n, mean = 0, sd = 1)

# linpred <- cbind(1, dummy(x)[, -1]) %*% c(beta0, betaB, betaC) + beta1 * xxx
linpred <- 
  dummy.mat[,'xA']*beta0 + 
  dummy.mat[,"xB"] * betaB + 
  dummy.mat[,"xC"] * betaC + 
  beta1 * xxx

prob_i <- exp(linpred) / (1 + exp(linpred))
y <- rbinom(n = n, size = 1, prob = prob_i)
data <- data.frame(x=x,xxx = xxx, y=y)
mod <- glm(y ~ x+xxx, family="binomial", data=data)
# summary(mod)

#------ initialisation ------
beta0Hat <- rep(NA, 1000)
betaBHat <- rep(NA, 1000)
betaCHat <- rep(NA, 1000)
#----------------------------

#------ simulations ------
# for(i in 1:100)
# {
#  #data generation
#  x <- sample(x=c("A","B", "C"), 
#              size=n, replace=TRUE, prob=rep(1/3, 3))  #(a)
#  linpred <- cbind(1, dummy(x)[, -1]) %*% c(beta0, betaB, betaC)  #(b)
#  prob_i <- exp(linpred) / (1 + exp(linpred))  #(c)
#  y <- rbinom(n=n, size=1, prob=prob_i)  #(d)
#  data <- data.frame(x=x, y=y)
#  
#  #fit the logistic model
#  mod <- glm(y ~ x, family="binomial", data=data)
#  
#  #save the estimates
#  beta0Hat[i] <- mod$coef[1]
#  betaBHat[i] <- mod$coef[2]
#  betaCHat[i] <- mod$coef[3]
# }
# #-------------------------
# 
# #------ results ------
# round(c(beta0 = mean(beta0Hat), 
#        betaB = mean(betaBHat), 
#        betaC = mean(betaCHat)), 3)



## ----glm_logistic_cat_plot, echo=FALSE, fig.width=6,fig.height=4----------------------------------------------------------------------------------------------
b0 <- mod$coef[1] # (Intercept)
xb <- mod$coef[2] # xB
xc <- mod$coef[3] # xC
XXX <- mod$coef[4] # xxx

xrange <- seq(from=min(data$xxx), to=max(data$xxx), length.out = n)

# Find the values for each category 
a_logits <- b0 + 
  XXX*xrange + 
  xb*0 + 
  xc*0 
b_logits <- b0 + 
  XXX*xrange + 
  xb*1 + 
  xc*0 
c_logits <- b0 + 
  XXX*xrange + 
  xb*0 + 
  xc*1 

# Transform
a_probs <- exp(a_logits)/(1 + exp(a_logits))
b_probs <- exp(b_logits)/(1 + exp(b_logits))
c_probs <- exp(c_logits)/(1 + exp(c_logits))

# plot 
plot(xrange, sort(a_probs), 
     ylim=c(0,1),
     type="l", lwd=3, lty=2, 
     col="gold", 
     xlab="x", ylab="P(outcome)", 
     main="Probability of outcome")


# Add the line for people who are in the b group
lines(xrange, sort(b_probs), 
      type="l", lwd=3, lty=3, 
      col="turquoise2")

# Add the line for people who are in the c group
lines(xrange, sort(c_probs), 
      type="l", lwd=3, lty=4, 
      col="orangered")
gold.5 = scales::alpha("gold", .3)
turquoise2.5 = scales::alpha("turquoise2", .3)
orangered.5 = scales::alpha("orangered", .3)
points(x = data[data$x %in% "A" & data$y %in% 1,"xxx"],data[data$x %in% "A" & data$y %in% 1,"y"], col = gold.5, pch =19, cex = .5)
points(x = data[data$x %in% "A" & data$y %in% 0,"xxx"],data[data$x %in% "A" & data$y %in% 0,"y"], col = gold.5, pch =19, cex = .5)

points(x = data[data$x %in% "B" & data$y %in% 1,"xxx"],data[data$x %in% "B" & data$y %in% 1,"y"]-.05, col = turquoise2.5, pch =19, cex = .5)
points(x = data[data$x %in% "B" & data$y %in% 0,"xxx"],data[data$x %in% "B" & data$y %in% 0,"y"]+.05, col = turquoise2.5, pch =19, cex = .5)

points(x = data[data$x %in% "C" & data$y %in% 1,"xxx"],data[data$x %in% "C" & data$y %in% 1,"y"]-.025, col = orangered.5, pch =19, cex = .5)
points(x = data[data$x %in% "C" & data$y %in% 0,"xxx"],data[data$x %in% "C" & data$y %in% 0,"y"]+.025, col = orangered.5, pch =19, cex = .5)

# add a horizontal line at p=.5
abline(h=.5, lty=2)


## ----sim_poisson_glm------------------------------------------------------------------------------------------------------------------------------------------
set.seed(42)
n = 500
x = rnorm(n = n, mean = 0, sd = 1)
# Rescale the data
xz = scale(x)
log.mu = 1 + .8*xz
y = rpois(n = n, lambda = exp(log.mu)) 

# Combine the data in a dataframe 
df = data.frame(y = y, x = x)


## ----sim_poisson_glm_plot, echo=-1, fig.width=9,fig.height=3.5------------------------------------------------------------------------------------------------
par(mfrow=c(1,3), mar=c(4,4,1.5,1.5), cex = 1.2)
#now feed it to glm:
glm.poisson = glm( y~x, data=df, family="poisson")
plot(y~x, data = df, col = scales::alpha("black", 0.5), pch = 19, cex = 0.5,
     main = "Poisson regression")
newdata <- data.frame(x = seq(min(x), max(x), len = n))
newdata$y = predict(object = glm.poisson, newdata = newdata, type = "response") 
lines(x = newdata$x,
      y = newdata$y, col = "red",lwd = 2)
glm.gau = glm( y~x, data=df, family="gaussian")
hist(residuals(glm.gau), main = "Residuals Gaussian")
hist(residuals(glm.poisson), main = "Residuals Poisson")


## ----sim_lmms, echo=FALSE, fig.width=5, fig.height=3----------------------------------------------------------------------------------------------------------
library(lme4)
set.seed(16)
# Experimental design 
lakes = 6
f.sp = 3
n.id = 10 
# Setting parameters
sd = 0.3# Individual measure 
sd.lake = .05
sd.fish = .02
beta0 = 1
beta1 = 0.003
total.n = lakes*f.sp*n.id

# Getting the age status of a fish 
age = rep(c(rep("Juv",2),rep("Ado",6),rep("Adt",2)), lakes)
n.juv = length(which(age =="Juv"))
n.ado = length(which(age =="Ado"))
n.adt = length(which(age =="Adt"))
age.df = data.frame(age,length = NA)
# Generating the length of the fish depending on the age it has 
age.df[age.df$age =="Juv","length"] <- rnorm(n.juv,mean = 100,sd = 10)
age.df[age.df$age =="Ado","length"] <- rnorm(n.ado,mean = 250,sd = 50)
age.df[age.df$age =="Adt","length"] <- rnorm(n.adt,mean = 400,sd = 10)

# trophic position is the response variable 
# Fish length is the phenotype that is measured 

lake.rep = gl(n = lakes,k = f.sp*n.id,labels = paste0("L",1:lakes))
# length(lake.rep)

sp.rep = rep( x = gl(n = f.sp,k = n.id, labels = paste0("s",1:f.sp)), lakes)
# length(sp.rep)
f.id = paste0("f.",1:(lakes*f.sp*n.id)) 

# Random effects 
# Setting it up with rnorm (the correct way)
# lake.rdm.eff = rep( rnorm(lakes, 0, sd.lake), each = f.sp*n.id)

# setting it up manually to see what happens when you modify the value
my.rdm.lake = c(1,2, 1.1,0.8,1.5,1.8)
lake.rdm.eff = rep( my.rdm.lake, each = f.sp*n.id)

# Setting it up with rnorm (the correct way)
# fish.rdm.eff = rep( rnorm(f.sp, 0, sd.fish), each = n.id)
# setting it up manually to see what happens when you modify the value
my.rdm.fish = c(-0.5,.4, -0.2)
fish.rdm.eff = rep( my.rdm.fish, each = n.id)

# Individual error 
id.err = rnorm(lakes*f.sp*n.id, 0, sd)

f.dat = data.frame(lake = lake.rep,
                   lake.rdm.eff,
                   species = sp.rep,
                   lake.rdm.eff,
                   fish.rdm.eff,
                   id = f.id,
                   age.df)

f.dat$t.lvl = with(f.dat, beta0 + beta1*length +lake.rdm.eff+fish.rdm.eff+ id.err )
# > range(fish.data$Trophic_Pos)
# [1] 2.123674 4.370899

f.dat$z_lgt = scale(f.dat$length) #(f.dat$Lake - mean(f.dat$length)) /sd(f.dat$length)
f.dat$z.t.lvl = scale(f.dat$t.lvl) #(f.dat$Lake - mean(f.dat$length)) /sd(f.dat$length)

# head(f.dat)
# range(f.dat$t.lvl)

plot <- ggplot(aes(length, t.lvl), data = f.dat)
fig <- theme_bw() + 
  theme(panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank(), 
        panel.background=element_blank(), 
        strip.background=element_blank(), 
        strip.text.y = element_text(),
        legend.background=element_blank(),
        legend.key=element_blank(),
        panel.border = element_rect(colour="black", fill = NA))

plot + geom_point(aes(col = species), size = 1) + 
  labs(x = "Length (mm)", y = "Trophic Position", 
       title = "All Data") + 
  fig

# Full model with varying intercepts and slopes only varying by species
M8 <- lmer(z.t.lvl ~ z_lgt + (1 + z_lgt | species) + (1 | lake),
           data = f.dat, REML = FALSE)

Lake.coef              <- coef(M8)$lake
colnames(Lake.coef)    <- c("Intercept", "Slope")
Species.coef           <- coef(M8)$species
colnames(Species.coef) <- c("Intercept", "Slope")


## ----sim_lmms_lines_species, echo=FALSE, fig.width=4, fig.height=3--------------------------------------------------------------------------------------------
plot <- ggplot(aes(z_lgt, z.t.lvl), data = f.dat)
Plot_BySpecies <- plot +
  geom_point(aes(colour = species), size = 1) +
  xlab("Length (mm)") + ylab("Trophic position") +
  labs(title = "By species") + fig
# Add regression lines for each species
Plot_BySpecies +
  geom_abline(intercept = Species.coef[1,1],
              slope     = Species.coef[1,2], col = "coral2") +
  geom_abline(intercept = Species.coef[2,1],
              slope     = Species.coef[2,2], col = "green4") +
  geom_abline(intercept = Species.coef[3,1],
              slope     = Species.coef[3,2], col = "blue1")



## ----sim_lmms_lines_lakes, echo=FALSE, fig.width=4, fig.height=3----------------------------------------------------------------------------------------------
Plot_ByLake <- plot +
  geom_point(aes(colour = lake), size = 1) +
  xlab("Length (mm)") + ylab("Trophic Position") +
  labs(title = "By Lake") + fig
# Add in regression lines with the intercepts specific to each lake
Plot_ByLake +
  geom_abline(intercept = Lake.coef[1,1],
              slope     = Lake.coef[1,2], col = "coral2") +
  geom_abline(intercept = Lake.coef[2,1],
              slope     = Lake.coef[2,2], col = "khaki4") +
  geom_abline(intercept = Lake.coef[3,1],
              slope     = Lake.coef[3,2], col = "green4") +
  geom_abline(intercept = Lake.coef[4,1],
              slope     = Lake.coef[4,2], col = "darkgoldenrod") +
  geom_abline(intercept = Lake.coef[5,1],
              slope     = Lake.coef[5,2], col = "royalblue1") +
  geom_abline(intercept = Lake.coef[6,1],
              slope     = Lake.coef[6,2], col = "magenta3")


## ----sim_lmms_data_gen----------------------------------------------------------------------------------------------------------------------------------------
set.seed(16)
# Experimental design 
lakes = 6
f.sp = 3
n.id = 10 

# Setting parameters
sd = 0.3# Individual measure 
sd.lake = .05
sd.fish = .02
beta0 = 1
beta1 = 0.003
total.n = lakes*f.sp*n.id

# Getting the age status of a fish 
age = rep(c(rep("Juv",2),rep("Ado",6),rep("Adt",2)), lakes)
n.juv = length(which(age =="Juv"))
n.ado = length(which(age =="Ado"))
n.adt = length(which(age =="Adt"))
age.df = data.frame(age,length = NA)


## ----sim_lmms_data_gen_randm----------------------------------------------------------------------------------------------------------------------------------
# Generating the length of the fish depending on the age it has 
age.df[age.df$age =="Juv","length"] <- rnorm(n.juv,mean = 100,sd = 10)
age.df[age.df$age =="Ado","length"] <- rnorm(n.ado,mean = 250,sd = 50)
age.df[age.df$age =="Adt","length"] <- rnorm(n.adt,mean = 400,sd = 10)

# trophic position is the response variable 
# Fish length is the phenotype that is measured 

lake.rep = gl(n = lakes,k = f.sp*n.id,labels = paste0("L",1:lakes))
# length(lake.rep)

sp.rep = rep( x = gl(n = f.sp,k = n.id, labels = paste0("s",1:f.sp)), lakes)
# length(sp.rep)
f.id = paste0("f.",1:(lakes*f.sp*n.id)) 

# Random effects 
# Setting it up with rnorm (the correct way)
# lake.rdm.eff = rep( rnorm(lakes, 0, sd.lake), each = f.sp*n.id)

# setting it up manually to see what happens when you modify the value
my.rdm.lake = c(1,2, 1.1,0.8,1.5,1.8)
lake.rdm.eff = rep( my.rdm.lake, each = f.sp*n.id)

# Setting it up with rnorm (the correct way)
# fish.rdm.eff = rep( rnorm(f.sp, 0, sd.fish), each = n.id)
# setting it up manually to see what happens when you modify the value
my.rdm.fish = c(-0.5,.4, -0.2)
fish.rdm.eff = rep( my.rdm.fish, each = n.id)

# Individual error 
id.err = rnorm(lakes*f.sp*n.id, 0, sd)

f.dat = data.frame(lake = lake.rep,
                   lake.rdm.eff,
                   species = sp.rep,
                   lake.rdm.eff,
                   fish.rdm.eff,
                   id = f.id,
                   age.df)

f.dat$t.lvl = with(f.dat, beta0 + beta1*length +lake.rdm.eff+fish.rdm.eff+ id.err )

f.dat$z_lgt = scale(f.dat$length) #(f.dat$Lake - mean(f.dat$length)) /sd(f.dat$length)
f.dat$z.t.lvl = scale(f.dat$t.lvl) #(f.dat$Lake - mean(f.dat$length)) /sd(f.dat$length)


## ----sim_lmms_data_gen_indv.err-------------------------------------------------------------------------------------------------------------------------------
# Individual error 
id.err = rnorm(lakes*f.sp*n.id, 0, sd)

f.dat = data.frame(lake = lake.rep,
                   lake.rdm.eff,
                   species = sp.rep,
                   lake.rdm.eff,
                   fish.rdm.eff,
                   id = f.id,
                   age.df)

f.dat$t.lvl = with(f.dat, beta0 + beta1*length +lake.rdm.eff+fish.rdm.eff+ id.err )

f.dat$z_lgt = scale(f.dat$length) #(f.dat$Lake - mean(f.dat$length)) /sd(f.dat$length)
f.dat$z.t.lvl = scale(f.dat$t.lvl) #(f.dat$Lake - mean(f.dat$length)) /sd(f.dat$length)


## ----lmm_sim.blog1, fig.width=8, fig.height=3-----------------------------------------------------------------------------------------------------------------
set.seed(16)
nstand = 5
nplot = 4
mu = 10
sds = 2
sd = 1
stand = rep(LETTERS[1:nstand], each = nplot) 
plot = letters[1:(nstand*nplot)] 
standeff = rnorm(nstand, 0, sds) 
standeff = rep(standeff, each = nplot) 
ploteff = rnorm(nstand*nplot, 0, sd) 
dat = data.frame(stand, standeff, plot, ploteff) 
dat$resp = with(dat, mu + standeff + ploteff ) 


## ----lmm_sim.blog2, fig.width=8, fig.height=3-----------------------------------------------------------------------------------------------------------------
library(lme4,warn.conflicts = FALSE)
fit1 = lmer(resp ~ 1 + (1|stand), data = dat)
fit1
twolevel_fun = function(nstand = 5, nplot = 4, mu = 10, sigma_s = 2, sigma = 1) {
     standeff = rep( rnorm(nstand, 0, sigma_s), each = nplot)
     stand = rep(LETTERS[1:nstand], each = nplot)
     ploteff = rnorm(nstand*nplot, 0, sigma)
     resp = mu + standeff + ploteff
     dat = data.frame(stand, resp)
     lmer(resp ~ 1 + (1|stand), data = dat)
}
set.seed(16)
twolevel_fun()
sims = replicate(100, twolevel_fun(), simplify = FALSE )
sims[[100]]
library(broom.mixed)
tidy(fit1)
tidy(fit1, effects = "fixed")
tidy(fit1, effects = "ran_pars", scales = "vcov")
library(purrr) # v. 0.3.4
suppressPackageStartupMessages( library(dplyr) ) # v. 1.0.7
library(ggplot2) # v. 3.3.5
stand_sims = c(5, 20, 100) %>%
     set_names() %>%
     map(~replicate(100, twolevel_fun(nstand = .x) ) )
stand_vars = stand_sims %>%
     modify_depth(2, ~tidy(.x, effects = "ran_pars", scales = "vcov") ) %>%
     map_dfr(bind_rows, .id = "stand_num") %>%
     filter(group == "stand")
head(stand_vars)


## ----lmm_sim.blog3, fig.width=8, fig.height=3-----------------------------------------------------------------------------------------------------------------
ggplot(stand_vars, aes(x = estimate) ) +
     geom_density(fill = "blue", alpha = .25) +
     facet_wrap(~stand_num) +
     geom_vline(xintercept = 4)
stand_vars = mutate(stand_vars, stand_num = forcats::fct_inorder(stand_num) )
add_prefix = function(string) {
     paste("Number stands:", string, sep = " ")
}
groupmed = stand_vars %>%
     group_by(stand_num) %>%
     summarise(mvar = median(estimate) )
ggplot(stand_vars, aes(x = estimate) ) + 
     geom_density(fill = "blue", alpha = .25) +
     facet_wrap(~stand_num, labeller = as_labeller(add_prefix) ) +
     geom_vline(aes(xintercept = 4, linetype = "True variance"), size = .5 ) +
     geom_vline(data = groupmed, aes(xintercept = mvar, linetype = "Median variance"),
                size = .5) +
     theme_bw(base_size = 14) +
     scale_linetype_manual(name = "", values = c(2, 1) ) +
     theme(legend.position = "bottom",
           legend.key.width = unit(.1, "cm") ) +
     labs(x = "Estimated Variance", y = NULL)
stand_vars %>%
     group_by(stand_num) %>%
     summarise_at("estimate", 
                  list(min = min, mean = mean, med = median, max = max) )
stand_vars %>%
     group_by(stand_num) %>%
     summarise(mean(estimate < 4) )


## ----lmmmmmmmmms, fig.width=8, fig.height=3-------------------------------------------------------------------------------------------------------------------
library(lmerTest, warn.conflicts = F)
n = 20; sd.n = 2
# Generate dataframe 
x = 1:n
values = rnorm(n = n,mean = 0,sd = sd.n)
gr = rep(c("short","tall"), each = n/2)
sim.df = data.frame(x,values,gr)

plot(density(sim.df[sim.df$gr%in%"short","values"]), col = "black", ylim = c(0,1), main = "Density")
lines(density(sim.df[sim.df$gr%in%"tall","values"]), col = "red")
legend("toprigh",legend = c("Short","Tall"),col = c("black","red"), lty = 1)


## ----help_arima, echo=FALSE, fig.width=8,fig.height=4---------------------------------------------------------------------------------------------------------
??arima.sim


## ----chaos_fun, echo=FALSE, fig.width=8,fig.height=4----------------------------------------------------------------------------------------------------------
chaos.fun <- function(lambda, n, initial, plot=TRUE, x.out = 20) {
  x <- numeric(n)
  x[1] <- initial
  for (t in 2 : n) x[t] <- lambda * x[t-1] * (1 - x[t-1])
  if(plot){ plot(1:n,x,type="l",ylim=c(0,1),ylab="population",
                 xlab="time",main=paste("lambda =", lambda))}
  tail(x,x.out)
}


## ----plot_chaos, echo=FALSE, fig.width=8,fig.height=4---------------------------------------------------------------------------------------------------------
par(mfrow=c(1,2))
chaout = chaos.fun(3.3,40,.6) # chaos.fun(4,40,.6)

plot(c(2,4),c(0,1),type="n",xlab="lambda",ylab="population")
for(lam in seq(2,4,0.001)){
  points(rep(lam,20), sapply(lam, FUN = function(x) {
    chaos.fun(lambda = x, n = 400, initial = .6, plot = FALSE)}),
    pch=16,cex=0.25,col=scales::alpha("blue", alpha = .5))
}



## ----plot_random_walk, fig.width=4,fig.height=4---------------------------------------------------------------------------------------------------------------
plot(0:100,0:100,type="n",xlab="",ylab="", main = "Random walk", asp = 1)
x <- y <- 50
points(50,50,pch=16,col="red",cex=1.5)
set.seed(1)
for (i in 1:1000) {
  xi <- sample(c(1,0,-1),1);   yi <- sample(c(1,0,-1),1) 
  lines(c(x,x+xi),c(y,y+yi),col=scales::alpha("blue", alpha = .7)) 
  x <- x+xi;   y <- y+yi
  if (x>100 | x<0 | y>100 | y<0) break }
points(x,y,pch=16,col="green",cex=1.5)


## ----rdm_points_in_polygon_code, eval=FALSE, fig.width=8, fig.height=3----------------------------------------------------------------------------------------
## library(sf); library(ggplot2)
## polygon = list(matrix(c(2, 2, 3, 3, 2.5, 4,
##                         3, 4, 2.5, 5, 1, 4,
##                         0, 5, 1, 3, 2, 2), ncol=2, byrow=T))
## polygon = sf::st_polygon(polygon) # Create an sf polygon
## points = sf::st_sample(polygon, size=50) # Sample 50 rdm pts in the polygon
## # Plot using the ggplot geom_sf function.
## pol.sf = ggplot() + geom_sf(aes(), data=polygon) +
##   geom_sf(aes(), col = alpha("black",.4), data=points) + theme_classic()


## ----rdm_points_in_polygon, echo=FALSE, fig.width=8, fig.height=3---------------------------------------------------------------------------------------------
library(sf); library(ggplot2)
polygon = list(matrix(c(2, 2, 3, 3, 2.5, 4, 
                        3, 4, 2.5, 5, 1, 4, 
                        0, 5, 1, 3, 2, 2), 
                      ncol=2, byrow=T)) 
polygon = sf::st_polygon(polygon) # Create an sf polygon
points = sf::st_sample(polygon, size = 50) # Sample 50 rdm pts in the polygon
# Plot using the ggplot geom_sf function.
pol.sf = ggplot() + geom_sf(aes(), data=polygon) + 
  geom_sf(aes(), col = alpha("black",.4), data=points) + theme_classic()

pol.sf + 
  theme(axis.text = element_text(colour = "black",size = 14),
    axis.text.x = element_text(colour = "black", size = 14),
    axis.text.y = element_text(colour = "black", size = 14),
    axis.line = element_line(size = 1.2,linetype = "solid"), 
    axis.ticks = element_line(colour = "black", size = 1.2), 
    panel.grid.major = element_line(colour = "gray98", linetype = "dashed"))


## ----spatial_read, echo=FALSE---------------------------------------------------------------------------------------------------------------------------------
library(sf, warn.conflicts = FALSE)
library(mapview, warn.conflicts = FALSE)

# Read the park layer
bot.gardp = st_read("data/GIS/Park_mtl.gpkg",layer = "Park_mtl", quiet = TRUE)
# Read the building layer
bot.gardb = st_read("data/GIS/Park_mtl.gpkg",layer = "buildings_mtl", quiet = TRUE)
# Show the layers 
# mapView(bot.gardb, col.regions = c("red")) + 
#   mapView(bot.gardp, col.regions = c("green"))

# Remove the buildings so that we can only have the 'green park' if we want to sample in the park 
only.park = st_difference(bot.gardp, st_union(bot.gardb))

# Get a specific point where we want to sample and get points around it  
selected.point = st_point(c(-73.566190,45.560516)) # Get point in CRS EPSG:4326
selected.point.no.tree = st_point(c(-73.56407,45.56502)) # Get point in CRS EPSG:4326
selected.point = st_sfc(selected.point) %>% st_set_crs(4326)
selected.point.no.tree = st_sfc(selected.point.no.tree) %>% st_set_crs(4326)

# Transform all the data to be in "NAD83 / MTM zone 8" or EPSG:32188
selected.point = st_transform(x = selected.point, crs = 32188)
selected.point.no.tree = st_transform(x = selected.point.no.tree, crs = 32188)
only.park = st_transform(x = only.park, crs = 32188)



## ----rdm_pts_in_polygon, echo=TRUE----------------------------------------------------------------------------------------------------------------------------
# Get random points
set.seed(456)
rdm.n = 5
rdm.pt = st_sample(x = only.park, 
                   size = rdm.n)
map.rdm = mapView(only.park, 
        col.regions = c("green")) +
  mapView(rdm.pt)  # Random points
# map.rdm

## save the output 
# mapshot(map.rdm,file = "images/map_rdm.png", 
        # url = "images/map_rdm.html")


## ----spatial_add_buffer---------------------------------------------------------------------------------------------------------------------------------------
# Add a buffer around the point we want to look at 
n = 10*2
min.buff = units::set_units(x = 0, value = "m")
max.buff = units::set_units(x = 100, value = "m")
buffer = st_buffer(x = c(selected.point,
                         selected.point.no.tree), 
                   dist = units::set_units(x = 100, 
                                           value = "m"))


## ----rdm_pts_in_polygon_buffer, echo=FALSE--------------------------------------------------------------------------------------------------------------------

map.unif = mapView(only.park, 
        col.regions = c("green"))+
  mapView(c(selected.point,selected.point.no.tree), # Get the poitn that was added to be looked at 
          col.regions = c("red")) +
  mapView(buffer, col.regions = c("red"))#; map.unif

## save html to png
# mapshot(map.unif,file = "images/map_unif.png", url = "images/map_unif.html")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(6543)
n = 10*2
min.buff = units::set_units(x = 0, value = "m")
max.buff = units::set_units(x = 100, value = "m")

# Get random distance 
# rdm.disttmp = rexp(n = n, rate = 1)*20
rdm.disttmp = runif(n = n, min = 0, max = max.buff)
# get random angle
rdm.angtmp = runif(n = n, min=0, max = 360)

# Conversion between Radians and Degrees
rad = rdm.angtmp * pi/180
rdm.ppt_x = rdm.disttmp*cos(rad) + c(st_coordinates(selected.point)[1], st_coordinates(selected.point.no.tree)[1])
rdm.ppt_y = rdm.disttmp*sin(rad) + c(st_coordinates(selected.point)[2], st_coordinates(selected.point.no.tree)[2])
rmd.ptdf = data.frame(rdm.ppt_x,
                      rdm.ppt_y, length(rdm.ppt_x))

rmd.ptdf.sf = sf::st_as_sf(rmd.ptdf, coords = c("rdm.ppt_x","rdm.ppt_y"), crs = 32188)#4326)


## ----rdm_pts_in_polygon_buffer_uniform, echo=FALSE------------------------------------------------------------------------------------------------------------
set.seed(456)
map.unif.rdm = mapView(only.park, col.regions = c("green"))+
  mapView(st_sample(only.park,5)) + 
  mapView(c(selected.point,selected.point.no.tree),col.regions = c("red")) + 
  mapView(buffer, col.regions = c("red")) + 
  mapview(rmd.ptdf.sf,col.regions = c("pink"))#; map.unif.rdm

## save html to png
# mapshot(map.unif.rdm,file = "images/map_unif_rdm.png", url = "images/map_unif_rdm.html")


## ----rdm_pts_in_polygon_grid, echo=TRUE-----------------------------------------------------------------------------------------------------------------------
set.seed(456)
# Add random points that are occupying the space of the polygon (grid )
rdm.pt = st_sample(x = only.park,
                   size = 100,
                   type ="hexagonal")
map.grid = mapView(only.park, 
                   col.regions = c("green")) + 
  mapView(rdm.pt)#; map.grid
## save html to png
# mapshot(map.grid,file = "images/map_grid.png", url = "images/map_grid.html")


## ----pwr_test-------------------------------------------------------------------------------------------------------------------------------------------------
p1 = .75 # Proportion to test 
p2 = .50 # proportion of the null hypothesis 
alpha = 0.05 # Type 1 error 
pwr = 0.80 # power or 1-beta = power. Beta is the type II error 
coin.p.power = power.prop.test(p1 = p1, p2 = p2, sig.level = alpha, power = pwr)
n = ceiling(coin.p.power$n) # get the sample size from the power analysis
coin.pip = rbinom(n, size = 1, prob = p1) # Generate coin tosses 
p.table = table(coin.pip)[c(2,1)] # Get the number of 1 and 0s  
(ptest = prop.test(p.table, alternative = "greater")) # Do the test 


## ----show_the_pwr, echo=-1, fig.width=7,fig.height=4----------------------------------------------------------------------------------------------------------
par(mar=c(4,4,0.1,0.1))
curve(dchisq(x, df = ptest$parameter), 
      xlim = c(0, ceiling(max(ptest$statistic))))
abline(v = ptest$statistic, lty = 3)


## ----get_eff_size, fig.width=7,fig.height=4-------------------------------------------------------------------------------------------------------------------
library(pwr)
r2 = seq(0,0.9,by =.1)
f2 <- function(r2) {r2/(1-r2)}


## ----get_eff_size_plot, echo=-1, fig.width=7,fig.height=4-----------------------------------------------------------------------------------------------------
par(mar=c(4,4,0.1,0.1))
plot(f2(r2)~r2)
curve(expr = (x/(1-x)),add = T)


## ----get_eff_size_n, fig.width=7,fig.height=4-----------------------------------------------------------------------------------------------------------------
nb.coef.in.model = 2
pwr.lm = pwr.f2.test(u = nb.coef.in.model, f2 = f2(.3), sig.level = 0.001, power = 0.8)

# sample size (n = v + u + 1)
n = ceiling(pwr.lm$v + nb.coef.in.model + 1)
n


## ----simulate_sampling_function-------------------------------------------------------------------------------------------------------------------------------
# Defining the population 
n = 600 # Number of elements to be generated 
set.seed(13) # Set RNG 
x = rnorm(n) # Generate numbers from normal distribution 
reedddd = scales::alpha("blue",.4) # Make colour transparent 

# Definte the function 
sample.mean.pop.est <- function(x,n, sample.size, ylim = NULL) {
  # x: is the actual values of the trait measured 
  # n: size of the population (number of individuals or items)
  # sample.size: how big is the sample size from which the MEAN will be calculated from 
  # ylim: add a maximum if needed 
  # histogram of the population 
  
  # Just get the stats from the histogram 
  pop.hist = hist(x, plot = F) # Make base histogram 
  
  # Make empty vector
  tmp.v = c(NA) 
  
  # For loop to calculate the mean based on a sample from the population 
  for (i in 1:n) {
    tmp = sample(x = x, size = sample.size, replace = FALSE)
    # Record that information (mean of the sample)
    tmp.v[i] = mean(tmp)
  } # End i
  
  # Sample histogram 
  sample.hist = hist(tmp.v, plot = F)
  # Population histogram 
  hist(x, ylim = range(c(0,c(sample.hist$counts, pop.hist$counts), ylim)), 
       main = paste("Sample n =", sample.size))
  # Add the sample estimate 
  sample.hist = hist(tmp.v, col = reedddd, add=T)
} # End sample.mean.pop.est


## ----simulate_sampling_plots----------------------------------------------------------------------------------------------------------------------------------
par(mfrow=c(2,2), lwd = .3)
sample.mean.pop.est(x = x, n = n, sample.size = 1, ylim = 300)
sample.mean.pop.est(x = x, n = n, sample.size = 10, ylim = 300)
sample.mean.pop.est(x = x, n = n, sample.size = 50, ylim = 300)
sample.mean.pop.est(x = x, n = n, sample.size = 500, ylim = 300)


## ----Dawkins_natural_Selection_example, echo=FALSE------------------------------------------------------------------------------------------------------------
set.seed(12345, kind="Mersenne-Twister")

## Easier if the string is a character vector
target <- unlist(strsplit("METHINKS IT IS LIKE A WEASEL", ""))
# target <- unlist(strsplit("MORE GIDDY IN MY DESIRES THAN A MONKEY", ""))
# http://shakespeare.mit.edu/hamlet/full.html 

charset <- c(LETTERS, " ")
rdm <- sample(charset, length(target), replace=TRUE)
parent <- sample(charset, length(target), replace=TRUE)

mutaterate <- 0.01

## Number of offspring in each generation
C <- 100

## Hamming distance between strings normalized by string length is used
## as the fitness function.
fitness <- function(parent, target) {
  sum(parent == target) / length(target)
}

mutate <- function(parent, rate, charset) {
  p <- runif(length(parent))
  nMutants <- sum(p < rate)
  if (nMutants) {
    parent[ p < rate ] <- sample(charset, nMutants, replace=TRUE)
  }
  parent
}

evolve <- function(parent, mutate, fitness, C, mutaterate, charset) {
  children <- replicate(C, mutate(parent, mutaterate, charset),
                        simplify=FALSE)
  children <- c(list(parent), children)
  children[[which.max(sapply(children, fitness, target=target))]]
}

.printGen <- function(parent, target, gen) {
  cat(format(i, width=3),
      formatC(fitness(parent, target), digits=2, format="f"),
      paste(parent, collapse=""), "\n")
}

i <- 0
.printGen(parent, target, i)
while ( ! all(parent == target)) {
  i <- i + 1
  parent <- evolve(parent, mutate, fitness, C, mutaterate, charset)
  
  if (i %% 20 == 0) {
    .printGen(parent, target, i)
  }
}
.printGen(parent, target, i)


## ----Dawkin_print, echo=FALSE---------------------------------------------------------------------------------------------------------------------------------
.printGen(rdm, target, i)

