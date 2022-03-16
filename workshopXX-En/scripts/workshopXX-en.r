## ----setup, echo = F------------------------------------------------------------------------------------------------------------------------------------------------------------------
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


## ----normal_compare_theoretical_simulated, echo=FALSE, fig.width=8, fig.height=4------------------------------------------------------------------------------------------------------
par(mfrow=c(1,2))
set.seed(12345)
curve(expr = dnorm(x), 
      from = -5,
      to = 5, 
      ylim = c(0,1), 
      xlab = "x", 
      ylab = "Density", 
      lwd =3)
plot(density(rnorm(6)), 
     main = "", 
     xlab = "x", 
     ylab = "Density", 
     xlim = c(-5,5), 
     ylim = c(0,1),
     lwd = 3)


## ----r_tips1, eval=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------
# # Description  ------------------------------------------------------------
# # This is the description section as previously presented
# 
# # Libraries ---------------------------------------------------------------
# # Here load the libraries used in the script
# library(ggplot2)
# 
# # Functions ---------------------------------------------------------------
# # Add and describe the functions used in the script
# ## Add more information
# function.name = function(var){tmp=2+var; return(tmp)}
# 
# # Plots -------------------------------------------------------------------
# # Plotting the data simulated
# ## Add more information
# plot(function.name(1:10)+rnorm(10))


## ----r_tips2, echo=FALSE, fig.width=13,fig.height=7-----------------------------------------------------------------------------------------------------------------------------------
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
par(mfrow = c(2,2), cex = 1.2)
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
       ylab = "Allele frq p",
       xlab = "Generations",
       ylim = c(0,1),
       xlim = c(1,gen),
       type = "n")
  
  # Add the lines per population and colour them 
  for (k in 1:popu) {
    pttmp = all.pops[all.pops$pop==k,]
    points(c(1:nrow(pttmp)), pttmp$p.fq, 
           type = "l",
           col = k)
  } # End k
} # End l



## ----fake_fitness_functions, echo=FALSE, fig.width=8,fig.height=3---------------------------------------------------------------------------------------------------------------------
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



## ----normalX_Y, echo=FALSE, fig.show='hide', fig.width=8,fig.height=4-----------------------------------------------------------------------------------------------------------------
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



## ----normalX_Y1, echo=FALSE, fig.width=5,fig.height=4---------------------------------------------------------------------------------------------------------------------------------
p.1


## ----normalX_Y2, echo=FALSE, fig.width=5,fig.height=4---------------------------------------------------------------------------------------------------------------------------------
p.2


## ----normalX_Y3, echo=FALSE, fig.width=5,fig.height=4---------------------------------------------------------------------------------------------------------------------------------
# p.h1
p.h.v1


## ----normalX_Y4, echo=FALSE, fig.width=5,fig.height=4---------------------------------------------------------------------------------------------------------------------------------
# p.h2
p.h.v2



## ----normalX_Y5, echo=FALSE, fig.width=5,fig.height=4---------------------------------------------------------------------------------------------------------------------------------
# p.h.v.sd.x1
p.h.v.sd.x.sd.y1


## ----normalX_Y6, echo=FALSE, fig.width=5,fig.height=4---------------------------------------------------------------------------------------------------------------------------------
# p.h.v.sd.x2
p.h.v.sd.x.sd.y2


## ----normalX_Y7, echo=FALSE, fig.width=5,fig.height=4---------------------------------------------------------------------------------------------------------------------------------
p.h.v.sd.x.sd.y.l1


## ----normalX_Y8, echo=FALSE, fig.width=5,fig.height=4---------------------------------------------------------------------------------------------------------------------------------
p.h.v.sd.x.sd.y.l2


## ----normalX_Y9, echo=FALSE, fig.width=5,fig.height=4---------------------------------------------------------------------------------------------------------------------------------
p.marg1


## ----normalX_Y10, echo=FALSE, fig.width=5,fig.height=4--------------------------------------------------------------------------------------------------------------------------------
p.marg2


## ----poissonX_Y, echo=FALSE, fig.show='hide', fig.width=8,fig.height=4----------------------------------------------------------------------------------------------------------------
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



## ----poissonX_Y1, echo=FALSE, fig.width=5,fig.height=4--------------------------------------------------------------------------------------------------------------------------------
p.1


## ----poissonX_Y2, echo=FALSE, fig.width=5,fig.height=4--------------------------------------------------------------------------------------------------------------------------------
p.2


## ----poissonX_Y3, echo=FALSE, fig.width=5,fig.height=4--------------------------------------------------------------------------------------------------------------------------------
p.marg1


## ----poissonX_Y4, echo=FALSE, fig.width=5,fig.height=4--------------------------------------------------------------------------------------------------------------------------------
p.marg2


## ----poissonX_logY, echo=FALSE, fig.show='hide', fig.width=8,fig.height=4-------------------------------------------------------------------------------------------------------------
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


## ----poissonX_logY3, echo=FALSE, fig.width=5,fig.height=4-----------------------------------------------------------------------------------------------------------------------------
p.marg1


## ----poissonX_logY4, echo=FALSE, fig.width=5,fig.height=4-----------------------------------------------------------------------------------------------------------------------------
p.marg2


## ----set.seed_function, echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(123)


## ----runif_example--------------------------------------------------------------------------------------------------------------------------------------------------------------------
runif(n = 1, min = 1, max = 10) # Gives a random number between 1 and 10
runif(n = 1, min = 1, max = 10) # RNG wasn't reset, different answer (see above)
runif(n = 1, min = 1, max = 10) # Different again... 

set.seed(42); runif(n = 1, min = 1, max = 10) # This sets the RNG 
set.seed(42); runif(n = 1, min = 1, max = 10) # The exact same number 



## ----set.seed_hidden, echo=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(123)


## ----sample_numerical_example---------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(12) # Set the RNG 
v.1.10 = 1:10 # Make a vector from 1 to 10 
# Randomly pick 1 (size) value from the vector (x), without replacement 
sample(x = v.1.10, size = 1, replace = FALSE) 


## ----sample_characters_example--------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(3) # Set the RNG 
# Randomly pick 5 (size) letters from the vector (x), without replacement 
sample(x = LETTERS, size = 5, replace = FALSE) 
sample(x = as.factor(month.abb), size = 5, replace = FALSE) 


## ----permutations_load_viridis, echo=FALSE--------------------------------------------------------------------------------------------------------------------------------------------
library(viridis)


## ----permutations_df, fig.width=4,fig.height=3----------------------------------------------------------------------------------------------------------------------------------------
set.seed(123)
n = 40; col = viridis::viridis(n = n)
x = 1:n ; y = 2+.5*x + rnorm(n, sd=7)
df.xy = data.frame(x,y, col )


## ----permutations_XY, fig.width=4,fig.height=3----------------------------------------------------------------------------------------------------------------------------------------
set.seed(321)
df.xy$x.s = sample(df.xy$x) #<<
df.xy$y.s = sample(df.xy$y) #<<
# We break up the link of X and Y 


## ----permutations_plot, fig.width=9,fig.height=3--------------------------------------------------------------------------------------------------------------------------------------
par(mfrow=c(1,3), mar=c(4,4,1,1), cex = 1.2)
plot(y~x,  col=col, data=df.xy, pch=19);abline(lm(y~x,  data=df.xy)) # Original 
plot(y~x.s,col=col, data=df.xy, pch=19);abline(lm(y~x.s,data=df.xy)) # Permutated 
plot(y.s~x.s,col=col, data=df.xy, pch=19);abline(lm(y~x.s,data=df.xy)) # Permutated 


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(42)
# GEt a sequence of dates
datae.seq = seq(from = as.Date('2010/01/01'), 
                to = as.Date('2022/01/01'), 
                by = "day")
# Look at beginning and end of sequence 
head(datae.seq, 4); tail(datae.seq, 4)
# Get only 5 elements of the generated sequence 
sample(datae.seq, 5)


## ---- echo=-1, fig.width=8, fig.height=3.5--------------------------------------------------------------------------------------------------------------------------------------------
par(mar=c(4,4,0.1,0.1))
set.seed(50)
p_dice = c(1,1,1,1,1,5) # Here we have twice the change of landing on 6
                        # Same as writing p_dice/sum(p_dice) or the prob.
nb.tosses = 100
die_results <- sample(x = 1:6, # or seq(from = 1, to=6, by=1)
                      size = nb.tosses,
                      replace = T, prob = p_dice) 
barplot(table(die_results)) # table(die_results)/nb.tosses


## ----rep_function_example-------------------------------------------------------------------------------------------------------------------------------------------------------------
(let4first = LETTERS[1:4])
rep(let4first, times = 2) # Repeat the WHOLE sequence twice 
rep(let4first, each = 2) # Repeat each element twice 
# Set the number of repeat for each element in the vector 
rep(let4first, times = c(1,2,3,6))
# Complete replication: replicate each element twice and do this three times 
rep(let4first, each = 2, times = 3)
rep(let4first, length.out = 6) # Repeat the vector until you hit the length.out



## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
nb.of.levels = 2
nb.of.replicates = 8
labels.for.the.factors = c("Control", "Treat")
## First control, then treatment:
gl(n = nb.of.levels, k = nb.of.replicates, labels = labels.for.the.factors)

## 20 alternating 1s and 2s
gl(n = 2, k = 1, length = 20)

## alternating pairs of 1s and 2s
gl(n = 2, k = 2, length = 20)


## ----data_replicate-------------------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(1234)
data.replicated = replicate(n = 2,
                            expr = data.frame(gr = rep(LETTERS[1:3], each = 2),
                                              y = rnorm(6)), 
                            simplify = FALSE)


## ----data_replicate1------------------------------------------------------------------------------------------------------------------------------------------------------------------
data.replicated[[1]]


## ----data_replicate2------------------------------------------------------------------------------------------------------------------------------------------------------------------
data.replicated[[2]]


## ----show_sample, eval=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------
# sample()


## ----sample_replacement---------------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(12345) # Sets the random number generator to a fix value
vec.1.10 = 1:10 # Make the vector to choose from 
sample(x = vec.1.10, size = 4, replace = FALSE) # Sample 4 nb without replacement
sample(x = vec.1.10, size = 4, replace = TRUE) # Sample 4 nb with replacement


## ----sample_example-------------------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(123); table(sample(x = 0:1, size = 1000, replace = T,prob = c(.5,.5)))


## ----linear_no_error, echo=-1, fig.width=6, fig.height=5------------------------------------------------------------------------------------------------------------------------------
par(mar=c(4,4,.5,.5), cex = 1.5)
x = 1:10
y = 2 + 3 * x
gr = rep(letters[1:2],each = 5)
linear.df = data.frame(x,y,gr)
plot(y~x, col = as.factor(gr), data = linear.df, pch = 19, cex = 1.5)


## ---- echo=FALSE, fig.width=4.7,fig.height=2------------------------------------------------------------------------------------------------------------------------------------------
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



## ---- echo=FALSE, fig.width=4.7,fig.height=2------------------------------------------------------------------------------------------------------------------------------------------
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


## ---- echo=FALSE, fig.width=4.7,fig.height=2------------------------------------------------------------------------------------------------------------------------------------------
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



## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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


## ----odds_example, echo=-1------------------------------------------------------------------------------------------------------------------------------------------------------------
library(MASS)
tot.nb.ev = 100; success = 20
odds.favor = (success)/(tot.nb.ev-success) # Odds in favor of event ( 1 in 4)
odds.agnst = (tot.nb.ev-success)/(success) # Odds against the event (here 4 to 1)


## ---- echo=FALSE, fig.width=4.7,fig.height=2------------------------------------------------------------------------------------------------------------------------------------------
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



## ----probability_example--------------------------------------------------------------------------------------------------------------------------------------------------------------
tot.nb.ev = 100; success = 20
probability.favor = (success)/(tot.nb.ev) # probability of event (here 20%)
probability.favor # fractions(probability.favor)


## ----odds_prob, fig.width=9,fig.height=7.5--------------------------------------------------------------------------------------------------------------------------------------------
# Get some LOG ODDS numbers 
log_odds = seq(from = -5, to = 5, by = .25)
# Transformed into odds 
odds = exp(log_odds)
# Make 
inv.logit <- function(x) {exp(x)/(1 + exp(x))}
p = inv.logit(log_odds) # This is the same as 
p2 = odds/(1 + odds)
# Probability of failure (1-p)
q = 1-p

# Store log_odds other data to plot 
d = data.frame(log_odds, odds, p, p2, q) 

head(signif(d,2), 3)


## ----odds_prob_plot, echo=-c(1), fig.width=12,fig.height=4----------------------------------------------------------------------------------------------------------------------------
par(mfrow = c(1,3), mar =c(4,4,.5,.5), cex = 1.3)
o.lo = d$odds~d$log_odds
p.o = d$p~d$odds
p.lo = d$p~d$log_odds
plot(o.lo, type="l", ylab="Odds", xlab="Ln odds", lwd=3); abline(v=0, lty=3)
plot(p.o,  type="l", ylab="p",    xlab="Odds",    lwd=3); abline(h=.5, v=0, lty=3)
plot(p.lo, type="l", ylab="p",    xlab="Ln odds", lwd=3); abline(h=.5, v=0, lty=3)


## ----feel_random_var, echo=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(123)
n = 1000
rnorm.val = rnorm(n = n, mean = 15, sd = 2)
random.var.normal = round(rnorm.val,1)


## ---- echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------
random.var.normal[1]


## ---- echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------
random.var.normal[1:5]


## ----feel_random_var_hist_5, echo=FALSE, fig.width=6, fig.height=5--------------------------------------------------------------------------------------------------------------------
par(mfrow=c(1,1),cex = 1.5)
hist(random.var.normal[1:5], main = "Hisogram of random variable", xlab = "x")


## ---- echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------
random.var.normal[1:100]


## ----feel_random_var_hist_100, echo=FALSE, fig.width=6, fig.height=5------------------------------------------------------------------------------------------------------------------
par(mfrow=c(1,1),cex = 1.5)
hist(random.var.normal[1:100], main = "Hisogram of random variable", xlab = "x")


## ----Statistical_dist, echo=FALSE, out.width=750, out.height=600, fig.width=9,fig.height=7.4------------------------------------------------------------------------------------------
par(mfrow = c(3,3), 
    # bg = NA,
    mar= c(4,4,3,3))
# par(mfrow = c(1,1), bg = NA)
col = scales::alpha("black",.5)
col.r = scales::alpha("red",.5)
col.g = scales::alpha("green",.5)
col.b = scales::alpha("blue",.5)
lwd=2
# Binom
#define range of "successes"
success <- 0:20
plot(success, dbinom(x = success, size=20, prob=0.4),type='h', col = col.b, lwd=lwd, ylim = c(0,1), ylab = "Density", main = "Binomial")
curve(expr = dbinom(x = x, size=20, prob=0.6),col = col, lwd=lwd, ylim = c(0,1), type = 'h',add = T)
curve(expr = dbinom(x = x, size=1, prob=0.5),col = col.r, lwd=lwd, ylim = c(0,1), type = 'h',add = T)
abline(h=seq(0,1, by = .1), lty = 3, lwd = .3)
legend("topright",legend = c("n = 20, p = .4","n = 20, p = .6","n = 1, p = .5"),lty = c(1,1,1), col =c(col.b,col, col.r), lwd = 2)

# Poisson
plot(success, dpois(success, lambda=5), type='h', col=col, lwd=lwd, ylim = c(0,1), ylab = "Density", main = "Poisson, lambda = 5")
abline(h=seq(0,1, by = .1), lty = 3, lwd = .3)

# Chi-sq
curve(dchisq(x, df = 10), from = 0, to = 40, col = col.b, lwd=lwd, ylim = c(0,1), ylab = "Density", main = "Chi-square")
curve(dchisq(x, df = 4), from = 0, to = 40, col = col, lwd=lwd, ylim = c(0,1), add = T)
curve(dchisq(x, df = 1), from = 0, to = 40, col = col.r, lwd=lwd, ylim = c(0,1), add = T)
abline(h=seq(0,1, by = .1), lty = 3, lwd = .3)
legend("topright",legend = c("df = 10","df = 4","df = 1"),lty = c(1,1,1), col =c(col.b,col, col.r), lwd = 2)

# Exponential 
curve(dexp(x, rate = .5), from=0, to=10, col=col.b, lwd=lwd, ylim = c(0,1), ylab = "Density", main = "Exponential")
curve(dexp(x, rate = .2), from=0, to=10, col=col, lwd=lwd, ylim = c(0,1),add = T)
curve(dexp(x, rate = .8), from=0, to=10, col=col.r, lwd=lwd, ylim = c(0,1),add = T)
abline(h=seq(0,1, by = .1), lty = 3, lwd = .3)
legend("topright",legend = c("rate = 0.8","rate = 0.5","rate = 0.2"),lty = c(1,1,1), col =c(col.b,col,col.r), lwd = 2)

# F-distribution
curve(df(x, df1 = 10, df2 = 20), from = 0, to = 4, n = 5000, col= col, lwd=lwd, ylim = c(0,1), ylab = "Density", main = "F-distribution, df1 = 10, df2 = 20")
abline(h=seq(0,1, by = .1), lty = 3, lwd = .3)

# normal 
curve(expr = dnorm(x = x, mean=0,sd=1), from = -5, to = 5, col=col.b, lwd=lwd, ylim = c(0,1), ylab = "Density", main = "Normal")
curve(expr = dnorm(x = x, mean=0,sd=2), from = -5, to = 5, col=col, lwd=lwd, ylim = c(0,1), add = T)
curve(expr = dnorm(x = x, mean=2,sd=1), from = -5, to = 5, col=col.r, lwd=lwd, ylim = c(0,1), add = T)
curve(dt(x, df=1), from=-5, to=5, col=col.g, lty = 1, lwd=lwd, ylim = c(0,1), ylab = "Density", main = "t-distribution df = 10", add=TRUE)
abline(h=seq(0,1, by = .1), lty = 3, lwd = .3)
legend("topright",legend = c("Normal, m = 0, sd = 1","Normal, m = 0, sd = 2","Normal, m = 2, sd = 1","t-distribution, df =1"),lty = c(1,1,1,1), col =c(col.b,col,col.r,col.g), lwd = 2)

# Log normal 
curve(dlnorm(x, meanlog=0, sdlog=1), from=0, to=10, col=col, lwd=lwd, ylim = c(0,1), ylab = "Density", main = "Log-normal, m = 0 sd = 1")
abline(h=seq(0,1, by = .1), lty = 3, lwd = .3)

# Logistic
curve(dlogis(x,location = 0, scale = 1), from=-10, to=10, col = col, lwd=lwd, ylim = c(0,1), ylab = "Density", main = "Logistic, l = 0, s = 1")
abline(h=seq(0,1, by = .1), lty = 3, lwd = .3)

# unifrom
curve(dunif(x, min = 8,max = 9), from=5, to=25, col=col.b, lwd=lwd, ylim = c(0,1), ylab = "Density", main = "Uniform") 
curve(dunif(x, min = 10,max = 15), from=5, to=25, col=col, lwd=lwd, ylim = c(0,1), add=T) 
curve(dunif(x, min = 16,max = 18), from=5, to=25, col=col.r, lwd=lwd, ylim = c(0,1), add = T) 
abline(h=seq(0,1, by = .1), lty = 3, lwd = .3)
legend("topright",legend = c("min = 8, max = 9","min = 10, max = 15","min = 16, max = 18"),lty = c(1,1,1), col =c(col.b,col,col.r), lwd = 2)


## ----similar_plot_from_different_distributions, echo=FALSE, fig.width=7,fig.height=3--------------------------------------------------------------------------------------------------
par(mfrow = c(1,2))
set.seed(1234)
hist(rbinom(10000, 10, 0.5), #breaks = seq(-0.5, 10.5, by = 1), 
     main = "Binomial")
hist(rnorm(10000, 5, 1.5), xlim = c(0,10), main = "Normal")


## ----equivalence_between_distributions_Bern_binom, echo=FALSE, fig.width=7,fig.height=4-----------------------------------------------------------------------------------------------
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


## ----simulate_coin_flips_plot, echo=FALSE, eval=FALSE---------------------------------------------------------------------------------------------------------------------------------
# # Coin flips visualize
# set.seed(1235)
# 
# nb.flips = 10
# set.seed(98765)
# coin.flips = rbinom(n = nb.flips, size = 1, prob = 0.5)
# heads.tails = ifelse(coin.flips==1,"H","T")
# radius =1
# # initialize a plot
# plot(x = c(0, nb.flips*(radius+.5)), y = c(-.5, 1), type = "n",
#      axes = F,
#      asp = 1,
#      ylab = "", xlab = "")
# w = 0
# pos=0
# for (i in 1:nb.flips) {
#   if (i %% 5 == 0) {
#     pos = pos + 2.5
#     w = 0
#   }
# 
#   # prepare "circle data"
#   radius = 1
#   center_x = w + 1
#   center_y = pos
#   theta = seq(0, 2 * pi, length = 200) # angles for drawing points around the circle
# 
#   # draw the circle
#   lines(x = radius * cos(theta) + center_x, y = radius * sin(theta) + center_y)
#   polygon(x = radius * cos(theta) + center_x, y = radius * sin(theta) + center_y, col = scales::alpha("beige",.5))
#   text(center_x,center_y,labels = heads.tails[i], cex = 3)
#   w = w + 2*radius+.5
# }
# 


## ----petri_dish_bacteria_colonny, echo=FALSE, fig.width=5,fig.height=3----------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(5937); rpois(n = 1,lambda = 2)


## ----Statistical_dist_poisson_wolf, echo=FALSE, fig.width=8,fig.height=2--------------------------------------------------------------------------------------------------------------
par(mfrow = c(1,1), mar = c(4,4,1,1))
set.seed(12345)
col = scales::alpha("black",.5); col.r = scales::alpha("red",.5)
col.g = scales::alpha("green",.5); col.b = scales::alpha("blue",.8)
lwd=2
add.l <- function(by=.1) {  abline(h=seq(0,1, by = by), lty = 3, lwd = .3)}
plot(0:10, dpois(x = 0:10, lambda=2),type='h', col = col.b, lwd=lwd, ylim = c(0,.41), ylab = "Probability", main = "lambda = 2"); add.l()


## ----Statistical_dist_poisson, echo=-c(1:6), fig.width=8,fig.height=4-----------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
dbinom(2,3,prob = .5) # 0.375


## ---- echo=-c(1:4), fig.width=8,fig.height=4------------------------------------------------------------------------------------------------------------------------------------------
par(mfrow = c(1,1), mar = c(4,4,1,1))
set.seed(12345); lwd = 3
col = scales::alpha("black",.5); col.r = scales::alpha("red",.5)
col.g = scales::alpha("green",.5); col.b = scales::alpha("blue",.8)
n = 7; x = 6; prob = 0.5
plot(0:8, dbinom(x = 0:8, size=n, prob=0.5),type='h', col = col.b, lwd=lwd, ylim = c(0,.4), ylab = "Probability", main = "s=7 and s=5; p=0.5 "); add.l(.05)
points(x = x, y = dbinom(x = x, size=n, prob=0.5), pch =19, col = col.b)

lines(c(0:8)+0.1, dbinom(x = 0:8, size=5, prob=0.5),type='h', col = col.r, lwd=lwd)
points(x = 3+.1, y = dbinom(x = 3, size=5, prob=0.5), pch =19, col = col.r)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(12345)
# 1000 experiments where each time, e.g., flipping a coin, I can either have a
table(rbinom(n = 1000, size=1, prob=0.5), dnn=NULL) # success or failure with p=.5
# 1 experiment where I have 1000 coins Where I sum all successes with p=.5
table(rbinom(n = 1, size=1000, prob=0.5), dnn=NULL) 
# 1000 experiments where each time, for example flipping 10 coins, where I 
table(rbinom(n = 1000, size=10, prob=0.5), dnn=NULL) # sum the success with p=.5


## ----Statistical_dist_binom, echo=-c(1:6), fig.width=8,fig.height=4-------------------------------------------------------------------------------------------------------------------
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


## ----equivalence_between_distributions_Binom_poisson, echo=FALSE----------------------------------------------------------------------------------------------------------------------
par(mfrow = c(2,2), mar = c(4,4,1,1))
x <- 0:10
n <- 10000
barplot(dbinom(x, n, 2/n), names.arg = x, ylim = c(0, 0.35), main = paste("Binomial with n = ", n))
barplot(dbinom(x, n, 9/n), names.arg = x, ylim = c(0, 0.35), main = paste("Binomial with n = ", n))
barplot(dpois(x, 2), names.arg = x, ylim = c(0, 0.35), main = paste("Poisson with Lambda = ", 2))
barplot(dpois(x, 9), names.arg = x, ylim = c(0, 0.35), main = paste("Poisson with Lambda = ", 9))

pbinom(q = 2, size = n, prob = 2/n)
ppois(q = 2, lambda = 2)


## ----Iris_normal, echo=FALSE, fig.width=12,fig.height=4-------------------------------------------------------------------------------------------------------------------------------
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


## ----Density_normal, fig.width=5,fig.height=5-----------------------------------------------------------------------------------------------------------------------------------------
curve(expr = dnorm(x = x, mean=0,sd=1), 
      from = -5, to = 5, ylim = c(0,1),
      col="black", ylab = "Density",
      main = "Density normal") 
abline(h=seq(0,1, by = .1), 
       lty = 3, lwd = .3)



## ----Density_normal_dissection_all_info, echo=FALSE, fig.width=5,fig.height=5---------------------------------------------------------------------------------------------------------
x <- seq(-5, 5, 0.1)
cex  = .7
plot(x, dnorm(x, 0, 1), main = "Density normal", type = "l", lwd = 3, col = "black", ylab = "", xlab = "x", ylim = c(0,1))
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
text(x = mean(x2), y = .2,labels = paste0((1-2*p)*100,"% \n=", "\npnorm(qnorm(1-p)) - \npnorm(qnorm(p))"),adj = 0,pos = 1, cex=cex)

text(x = 0, y = .55,labels = paste0(quote("Mean")),adj = 0,pos = 1,offset = 0, cex=cex)
arrows(x0 = 0,x1 = 0, y0 = .5,y1 = .42,code = 2,length=.1)

legend("topright",legend = c("Density normal (dnorm)","Density rnorm(100)"), lwd = 1, lty = 1, col = c("black", "green"))

set.seed(123)
rndat =rnorm(100)
# mean(rndat)
dens.nor = density(rndat)
lines(dens.nor, lwd = 3,col = scales::alpha("green",.8))
mybins=hist(rndat, plot = F, density = T, breaks = 100)
crn = mybins$density
brn = mybins$breaks


## ----pnorm_qnorm2---------------------------------------------------------------------------------------------------------------------------------------------------------------------
pnorm(1.645)

qnorm(p = 0.05, lower.tail = F)

qnorm(p = 0.025, lower.tail = F)


## ----normal_shade_95_5, echo=FALSE, fig.width=5,fig.height=5--------------------------------------------------------------------------------------------------------------------------
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


## ----normal_area_function, echo=FALSE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------------------
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


## ----normal_dist_area, echo=-1, fig.width=9,fig.height=5------------------------------------------------------------------------------------------------------------------------------
par(mfrow=c(2,3), mar = c(4,4,.1,.1), cex = 1.1)
# The function is not shown, but can be found in the script (Markdown)
draw.normal(where = "both",  prob = 0.05/2)
draw.normal(where = "both",  prob = 0.2/2 )
draw.normal(where = "both",  prob = 0.5/2 )
draw.normal(where = "both",  prob = 0.95/2)
draw.normal(where = "left",  prob = 0.05  )
draw.normal(where = "right", prob = 0.05  )


## ----Normal_pdf_important_values------------------------------------------------------------------------------------------------------------------------------------------------------
sd = 1
probability.left.side = (pnorm(q = c(sd*1,sd*2,sd*3),lower.tail = F)*100)
probability.right.side = (pnorm(q = c(sd*1,sd*2,sd*3),lower.tail = T)*100)
percent.data.under.curve = probability.right.side - probability.left.side
p.from.mean = round(percent.data.under.curve,2)


## ----normal_dist_area2, echo=FALSE, fig.width=11,fig.height=3-------------------------------------------------------------------------------------------------------------------------
par(mfrow=c(1,3), mar = c(4,4,1,1), cex = 1.1)
draw.normal(where = "both",  prob = round(probability.left.side[1]/100,3))
draw.normal(where = "both",  prob = round(probability.left.side[2]/100,3))
draw.normal(where = "both",  prob = round(probability.left.side[3]/100,3))


## ----Normal_pdf_important_values2-----------------------------------------------------------------------------------------------------------------------------------------------------
qnorm(p = c(.75, .95,.975, .995), mean = 0, sd = 1, lower.tail = F)
qnorm(p = c(.75, .95,.975, .995), mean = 0, sd = 1, lower.tail = T)


## ----normal_dist_area3, echo=FALSE, fig.width=14,fig.height=3-------------------------------------------------------------------------------------------------------------------------
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



## ---- echo=FALSE, fig.width=7,fig.height=5--------------------------------------------------------------------------------------------------------------------------------------------
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


## ---- echo=FALSE, fig.width=7,fig.height=5--------------------------------------------------------------------------------------------------------------------------------------------
par(mfrow=c(1,1), mar = c(4,4,0,1), cex = 1.1)
hypothesis.testing(mean.pop = 2)


## ---- echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------
par(mfrow=c(1,1), mar = c(4,4,0,1), cex = 1.1)
hypothesis.testing(2,sd = 5, n = 10)


## ---- echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------
par(mfrow=c(1,1), mar = c(4,4,0,1), cex = 1.1)
hypothesis.testing(2,sd = 5, n = 100)


## ----Standard_normal_transformation, echo=FALSE, fig.width=10,fig.height=5------------------------------------------------------------------------------------------------------------
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



## ----rnom_function, eval=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------
# rnorm()


## ----rnom_function_example------------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(1234)
n <-10
rnorm(n)


## ---- echo=-1, fig.width=8,fig.height=4.5---------------------------------------------------------------------------------------------------------------------------------------------
par(mfrow=c(2,2), cex =1.1, mar = c(4,4,1,1)) # set window 
n = 1000 # Number of points 
# Generate multiple additions of random variables 
for(i in c(2, 50, 1000, 5000)){
  clt = replicate(i, rexp(n, rate = 1), simplify = FALSE)
  hist(apply(do.call(cbind,clt),1,sum), main = paste("Hist. of",i,"variables"), xlab = "x") # Draw the histogram 
}


## ---- echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Define uniform discrete 
dunifdisc<-function(x, min=0, max=1) ifelse(x>=min & x<=max & round(x)==x, 1/(max-min+1), 0)
punifdisc<-function(q, min=0, max=1) ifelse(q<min, 0, ifelse(q>=max, 1, (floor(q)-min+1)/(max-min+1)))
qunifdisc<-function(p, min=0, max=1) floor(p*(max-min+1))
runifdisc<-function(n, min=0, max=1) sample(min:max, n, replace=T)


## ----CLT_uniform_example_dice, echo=-1, eval=FALSE, fig.width=8,fig.height=5----------------------------------------------------------------------------------------------------------
# par(mfrow=c(2,2)) # set window
# # Generate multiple additions of random variables
# for(i in c(2, 50, 1000, 5000)){
#   clt = replicate(i, runifdisc(n),simplify = FALSE)
#   hist(apply(do.call(cbind,clt),1,sum), main = paste("Histogram of",i,"variables"),xlab = "x") # Draw the histogram
# }


## ----CLT_uniform_example_dice_hist, echo=FALSE, fig.width=8,fig.height=5--------------------------------------------------------------------------------------------------------------
par(mfrow=c(2,2)) # set window 
# Generate multiple additions of random variables 
for(i in c(2, 50, 1000, 5000)){
  clt = replicate(i, runifdisc(n),simplify = FALSE)
  hist(apply(do.call(cbind,clt),1,sum), main = paste("Histogram of",i,"variables"),xlab = "x") # Draw the histogram 
}


## ---- echo=FALSE, fig.width=3,fig.height=3--------------------------------------------------------------------------------------------------------------------------------------------
curve(dunifdisc(x, 7,10), type = "h",from=6, to=11, col="black", lwd=1, ylim = c(0,1), ylab = "Density", main = "Uniform") 


## ---- echo=-1, fig.width=8,fig.height=3-----------------------------------------------------------------------------------------------------------------------------------------------
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


## ---- echo=FALSE, fig.width=8,fig.height=4--------------------------------------------------------------------------------------------------------------------------------------------
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



## ---- echo=-1, fig.width=5,fig.height=2-----------------------------------------------------------------------------------------------------------------------------------------------
par(mfrow = c(1,2), mar = c(4,4,0.5,0.5)); b.5 = scales::alpha("black",0.5)
set.seed(24601) # setting this so the random results will be repeatable 
library(MASS)
covmat <- matrix(c(1.0,   0.2,   0.6, # variance covariance matrix of the data
                   0.2,   2.0,  -0.5, 
                   0.6,  -0.5,   1.0), nrow=3) 
data <- mvrnorm(n = 300, 
                mu = c(1,-1,0), # mean of the data 
                Sigma=covmat) # generate random data that match that variance covariance matrix
plot(data[,1:2], pch = 19, col = b.5); abline(h=0,v=0,lty = 3)
plot(data[,2:3], pch = 19, col = b.5); abline(h=0,v=0,lty = 3)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
x = rnorm(100)
# x = rt(100,1000)
library(fitdistrplus)
descdist(x, discrete = FALSE)
fit.pois <- fitdist(x, "norm")
plot(fit.pois)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
p1 = .75 # Proportion to test 
p2 = .50 # proportion of the null hypothesis 
alpha = 0.05 # Type 1 error 
pwr = 0.80 # power or 1-beta = power. Beta is the type II error 
coin.p.power = power.prop.test(p1 = p1, p2 = p2, sig.level = alpha, power = pwr)
n = ceiling(coin.p.power$n) # get the sample size from the power analysis
coin.pip = rbinom(n, size = 1, prob = p1) # Generate coin tosses 
p.table = table(coin.pip)[c(2,1)] # Get the number of 1 and 0s  
(ptest = prop.test(p.table, alternative = "greater")) # Do the test 


## ---- echo=-1, fig.width=7,fig.height=4-----------------------------------------------------------------------------------------------------------------------------------------------
par(mar=c(4,4,0.1,0.1))
curve(dchisq(x, df = ptest$parameter), 
      xlim = c(0, ceiling(max(ptest$statistic))))
abline(v = ptest$statistic, lty = 3)


## ---- fig.width=7,fig.height=4--------------------------------------------------------------------------------------------------------------------------------------------------------
library(pwr)
r2 = seq(0,0.9,by =.1)
f2 <- function(r2) {r2/(1-r2)}


## ---- echo=-1, fig.width=7,fig.height=4-----------------------------------------------------------------------------------------------------------------------------------------------
par(mar=c(4,4,0.1,0.1))
plot(f2(r2)~r2)
curve(expr = (x/(1-x)),add = T)


## ---- fig.width=7,fig.height=4--------------------------------------------------------------------------------------------------------------------------------------------------------
nb.coef.in.model = 2
pwr.lm = pwr.f2.test(u = nb.coef.in.model, f2 = f2(.3), sig.level = 0.001, power = 0.8)

# sample size (n = v + u + 1)
n = ceiling(pwr.lm$v + nb.coef.in.model + 1)
n


## ----Sim_t_test_anova-----------------------------------------------------------------------------------------------------------------------------------------------------------------
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


## ----Sim_t_test_anova_plot, fig.width=5,fig.height=5----------------------------------------------------------------------------------------------------------------------------------
plot(y~gr, data = df.aov)


## ---- echo=-1, fig.width=4,fig.height=2-----------------------------------------------------------------------------------------------------------------------------------------------
par(mar=c(4,4,0.1,0.1))
set.seed(12345678)
n = 100; beta0 = 2.5; beta1 = 0.8
x.lm = rnorm(n = n, mean = 10, sd = 1)
err = rnorm(n = n, mean = 0, sd = 1)
# Linear combination 
y.lm = beta0 + beta1*x.lm + err
# Make a dataframe of the data 
df.lm = data.frame(x = x.lm, y = y.lm)


## ---- echo=FALSE, fig.width=4,fig.height=3--------------------------------------------------------------------------------------------------------------------------------------------
par(mar = c(4,4,.5,.5))
# Colour 
b.5 = scales::alpha("black",alpha = .5)

# PLot the data 
plot(y~x, data = df.lm, pch = 19, col = b.5)

# Model the data 
lm.out = lm(y~x, data = df.lm)
# Add a line to the plot 
abline(lm.out)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
round(coefficients(summary(lm.out)), 4)
# Adjusted R^2
summary(lm.out)$adj.r.squared
# summary(lm.out)$fstatistic
# anova(lm.out)


## ---- echo=-1, fig.width=8,fig.height=3-----------------------------------------------------------------------------------------------------------------------------------------------
par(mfrow = c(1,3), mar = c(4,4,1,1), cex = 1.2)
sim.lm = simulate(lm.out, nsim = 2000, seed = 12)
r.x = range(c(y.lm, rowMeans(sim.lm), fitted(lm.out)))
hist(rowMeans(sim.lm), xlim = r.x, main = "Hist simulation")
hist(fitted(lm.out), xlim = r.x, main = "Hist fitted")
hist(y.lm, xlim = r.x, main = "Hist of response")
c(mean(rowMeans(sim.lm)), mean(x.lm),  mean(fitted(lm.out))) # compare 
rbind(head(rowMeans(sim.lm)), head(fitted(lm.out)))


## ---- echo=FALSE, fig.width=8,fig.height=3--------------------------------------------------------------------------------------------------------------------------------------------
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
  
  # Return some interesting parameters for futur analysis
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
                     data.frame(n = i, all.int, all.slp))
}

plot.dens <- function(density, data, subset, main) {
  plot(density,   xlim = c(-3,4), main = main)
abline(v = mean(data[data$n==subset,"all.slp"]),   lty =1)
abline(v = beta1, lty = 3)
polygon(x = c(-10, density$x[density$x>-10 & density$x < 10], 10), 
        y = c(0, density$y[density$x>=-10 & density$x <= 10], 0),
        col=scales::alpha("blue",.5))
}

par(mfrow = c(1,3))
# get the density 
dens5 = density(df.all.sim[df.all.sim$n==5,"all.slp"]  )
dens20 = density(df.all.sim[df.all.sim$n==20,"all.slp"] )
dens200 = density(df.all.sim[df.all.sim$n==200,"all.slp"])

# Plot for each parameter we are interested in 
plot.dens(dens5,df.all.sim,5,"n = 5")
plot.dens(dens20,df.all.sim,5,"n = 5")
plot.dens(dens200,df.all.sim,5,"n = 5")

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


## ---- echo=FALSE, fig.width=10,fig.height=5-------------------------------------------------------------------------------------------------------------------------------------------
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



## ---- echo=FALSE, fig.width=8,fig.height=4--------------------------------------------------------------------------------------------------------------------------------------------
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


## ---- echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Coefficients table:\n"); round(coefficients(summary(lm.out)), 4)
# Adjusted R^2
cat("R Squared: "); summary(lm.out)$adj.r.squared
# summary(lm.out)$fstatistic
# anova(lm.out)
# vcov(lm.out)


## ----sim_lm_cat, fig.width=5,fig.height=5---------------------------------------------------------------------------------------------------------------------------------------------
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


## ---- echo=FALSE, eval=TRUE, results = 'hide'-----------------------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Here you can see the beta1 (x), betaControl (Intercept) and betaTreament (grTreatment)
summary(lm(y~x+gr, data = df.lm))
# Get JUST the DIFFERENCE with the treatment (Intercept is the Control)
# gr.means = lm(y~gr, data = df.lm) 
# coefficients(gr.means)


## ----sim_lm_cat_plot, echo=-1, fig.width=8,fig.height=3.5-----------------------------------------------------------------------------------------------------------------------------
par(mfrow=c(1,2))
plot(y~x, data = df.lm, pch = 19, col = ifelse(gr=="Treatment",r.5,b.5))
abline(a= lmcoef["(Intercept)"], b= lmcoef["x"], col = "black") # control
abline(a= lmcoef["(Intercept)"] + lmcoef["grTreatment"], b=lmcoef["x"], col="red")
boxplot(y~gr, data = df.lm, pch = 19, col = ifelse(gr=="Treatment",r.5,b.5))
round(coefficients(summary(lm.out)), 4) # betaControl is (Intercept)
summary(lm.out)$adj.r.squared           # betaTreament is grTreatment. 


## ---- echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------
df.eo


## ----sim_lm_interaction_discrete, echo=FALSE, fig.width=8,fig.height=3----------------------------------------------------------------------------------------------------------------
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


## ---- echo=-1, fig.width=5,fig.height=5-----------------------------------------------------------------------------------------------------------------------------------------------
par(mfrow=c(1,1), mar = c(4,4,.1,.1))
plot(x = df.lm[df.lm$x1 == 0, ]$x2, y = df.lm[df.lm$x1 == 0, ]$y, # Add x1 = 0
     col = rgb(red = 0, green = 0, blue = 1, alpha = 0.25), pch = 19,
     xlab = "x2", ylab = "y", ylim = range(df.lm$y)); abline(h=0, v=0, lt =3)
abline(a = coef(lm.out)[1], b = coef(lm.out)[3], col = "blue", pch = 19, lwd =2)
points(x = df.lm[df.lm$x1 == 1, ]$x2, y = df.lm[df.lm$x1 == 1, ]$y, # Add x1 = 1
       col = rgb(red = 1, green = 0, blue = 0, alpha = 0.25), pch = 19)
abline(a = coef(lm.out)[1] + coef(lm.out)[2], 
       b = coef(lm.out)[3] + coef(lm.out)[4], col = "red", lwd = 2)


## ----sim_lm_interaction, echo=FALSE, fig.width=8,fig.height=5-------------------------------------------------------------------------------------------------------------------------
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


## ---- echo=-c(1:2), fig.width=8,fig.height=2,.5---------------------------------------------------------------------------------------------------------------------------------------
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


## ---- echo=-c(1:3), fig.width=8,fig.height=4------------------------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
summary(lm.out)


## ----challenge4.one.rep, fig.width=5,fig.height=5-------------------------------------------------------------------------------------------------------------------------------------
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


## ----challenge4.diagnose1, echo=-c(1), fig.width=5,fig.height=4-----------------------------------------------------------------------------------------------------------------------
par(mfrow = c(1,1), mar = c(4, 4, .1, 0.1))
plot(weight~gr, data=plant.weight, las=1)


## ----challenge4.diagnose2, echo=-c(1), fig.width=5,fig.height=5-----------------------------------------------------------------------------------------------------------------------
par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
plot(lm.out.int, las = 1)


## ----challenge4.fun, fig.width=8,fig.height=5-----------------------------------------------------------------------------------------------------------------------------------------
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



## ----challenge4.fun.plot, echo=-c(1:3), fig.width=8,fig.height=3.5--------------------------------------------------------------------------------------------------------------------
par(mfrow = c(2,2), mar = c(4,4,.5,.5))
layout.matrix <- matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2)
layout(layout.matrix, widths = c(2,2), heights = c(1,1))
set.seed(2345); nb.rep = 2000 # number of replications in the simulation 
l.rp = replicate(n = nb.rep, simplify = FALSE,
                 expr = exp.plant( n = 10)$coefficients["grtrt","t value"])
p.val.lm = pt(q = abs(unlist(l.rp)),df = s.lm.out$df[2], lower.tail = F)*2 # Get p-values 
# plots 
exp.plant(n = 10, plot = T, ret = F) # plot a simulation example 
hist(unlist(l.rp), main = "t-values", probability = T, xlab ="")
lines(density(unlist(l.rp)), col="blue", lwd=2); abline(v = qt(0.025, df = s.lm.out$df[2]))
hist(p.val.lm, main = "p-values", probability = T, xlab ="", xlim = c(0,1))
lines(density(p.val.lm), col="red", lwd=2); abline(v = 0.05)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
qt(0.05, df = s.lm.out$df[2])
sum(unlist(l.rp)<qt(0.025, df = s.lm.out$df[2]))/length(unlist(l.rp))
sum(p.val.lm<=0.05)/length(p.val.lm)


## ----sim_logistic---------------------------------------------------------------------------------------------------------------------------------------------------------------------
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



## ----sim_logistic_glm, echo=-1, fig.width=5,fig.height=3------------------------------------------------------------------------------------------------------------------------------
par(mar = c(4,4,0.2,0.1))
plot(y~x1, data = df, col = scales::alpha("black",.5), pch = 19)
newdata <- data.frame(x1=seq(min(x1), max(x1),len=n), 
                      x2 = seq(min(x2), max(x2),len=n))
newdata$y = predict(object = glm.logist, newdata = newdata, type = "response") 
lines(x = newdata$x1, y = newdata$y, col = "red",lwd = 2)


## ---- echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------
glm.sum = summary(glm.logist) # Note that the output of the GLM is in "logit"
coefficients(glm.sum)[,"Estimate"]
# glm.sum$coefficients[,1]/glm.sum$coefficients[,2] # Estimate/Std. Error = z val.


## ---- echo=FALSE, fig.width=8,fig.height=3--------------------------------------------------------------------------------------------------------------------------------------------
set.seed(6)
n = 20
x1 = rnorm(n = n, mean = 6, sd = 1)
# Rescale the data
x1z = scale(x1)
z = 0 + 2*x1z  
pr = 1/(1+exp(-z))
y = rbinom(n = n, size = 1, prob = pr) 

# Combine the data in a dataframe 
df = data.frame(y = y, x1 = x1)

#now feed it to glm:
glm.logist = glm( y~x1, data=df, family="binomial")
glm.sum = summary(glm.logist)
# cat("Coefficients\n")
glm.sum$coefficients

par(mfrow=c(1,3))
b.5 = scales::alpha("black",.5)
plot(z~x1z, ylab = "Log Odds", main = "Theoretical logistic in log(odds)",
     pch = 19, col = b.5, xlim = c(-5,5), ylim = c(-12,11))
abline(a = 0,
       b = 2, col = "red")
abline(h=0, v=0,lty = 3)

plot(z~x1, ylab = "Log Odds", main = "Estimated logistic in log(odds)",
     pch = 19, col = b.5, xlim = c(0,10), ylim = c(-12,11))
abline(a = glm.sum$coefficients[1,1],
       b = glm.sum$coefficients[2,1])
abline(h=0, v=0,lty = 3)
points(x = 0, y=glm.sum$coefficients[1,1], pch = 19, col = "red")
text(x = 0, y=glm.sum$coefficients[1,1], labels = c("Intercept"), pos =4)

plot(y~x1, data = df, main = "Estimated logistic in probability",
     ylab = "Probability of outcome",
     col = scales::alpha("black",.5), pch = 19)
abline(h=0.5, v=mean(x1),lty = 3)
newdata <- data.frame(x1=seq(min(x1), max(x1),len=n))
newdata$y = predict(object = glm.logist, newdata = newdata, type = "response") 
lines(x = newdata$x1,
      y = newdata$y, col = "red",lwd = 2)

cat("Probability at intercept is",1/(1+exp(-glm.sum$coefficients[1,1])),"\n")
cat("Probability for an increase in 1 sd of X is",1/(1+exp(-glm.sum$coefficients[2,1])))


## ---- echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------
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
beta0 <- 0.07
betaB <- 0.1
betaC <- -0.15
beta1 <- .5

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



## ---- echo=FALSE, fig.width=6,fig.height=4--------------------------------------------------------------------------------------------------------------------------------------------
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

# add a horizontal line at p=.5
abline(h=.5, lty=2)


## ----sim_poisson_glm------------------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(42)
n = 1000
x = rnorm(n = n, mean = 0, sd = 1)
# Rescale the data
xz = scale(x)
log.mu = 1 + 2*xz
y = rpois(n = n, lambda = exp(log.mu)) 

# Combine the data in a dataframe 
df = data.frame(y = y, x = x)


## ----sim_poisson_glm_plot, fig.width=5,fig.height=5-----------------------------------------------------------------------------------------------------------------------------------
#now feed it to glm:
glm.poisson = glm( y~x, data=df, family="poisson")
plot(y~x, data = df, col = scales::alpha("black",.5), pch = 19)
newdata <- data.frame(x=seq(min(x), max(x),len=n))
newdata$y = predict(object = glm.poisson, newdata = newdata, type = "response") 
lines(x = newdata$x,
      y = newdata$y, col = "red",lwd = 2)


## ---- echo=FALSE, fig.width=5, fig.height=3-------------------------------------------------------------------------------------------------------------------------------------------
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


## ---- echo=FALSE, fig.width=4, fig.height=3-------------------------------------------------------------------------------------------------------------------------------------------
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



## ---- echo=FALSE, fig.width=4, fig.height=3-------------------------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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


## ---- fig.width=8, fig.height=3-------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(16)
nstand = 5
nplot = 4
mu = 10
sds = 2
sd = 1
( stand = rep(LETTERS[1:nstand], each = nplot) )
( plot = letters[1:(nstand*nplot)] )
( standeff = rnorm(nstand, 0, sds) )
( standeff = rep(standeff, each = nplot) )
( ploteff = rnorm(nstand*nplot, 0, sd) )
( dat = data.frame(stand, standeff, plot, ploteff) )
( dat$resp = with(dat, mu + standeff + ploteff ) )
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



## ----lmmmmmmmmms, fig.width=8, fig.height=3-------------------------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
??arima.sim


## ---- fig.width=8, fig.height=3-------------------------------------------------------------------------------------------------------------------------------------------------------
library(sf); library(ggplot2)
polygon = list(matrix(c(2, 2, 3, 3, 2.5, 4, 
                        3, 4, 2.5, 5, 1, 4, 
                        0, 5, 1, 3, 2, 2), ncol=2, byrow=T)) 
polygon = sf::st_polygon(polygon) # Create an sf polygon
points = sf::st_sample(polygon, size=50) # Sample 50 rdm pts in the polygon
# Plot using the ggplot geom_sf function.
ggplot() + geom_sf(aes(), data=polygon) + 
  geom_sf(aes(), col = alpha("black",.4), data=points) + theme_classic()


## ----echo=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
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
selected.point = st_sfc(selected.point) %>% 
  st_set_crs(4326)
selected.point.no.tree = st_sfc(selected.point.no.tree) %>% 
  st_set_crs(4326)

# Transform all the data to be in "NAD83 / MTM zone 8" or EPSG:32188
selected.point = st_transform(x = selected.point, crs = 32188)
selected.point.no.tree = st_transform(x = selected.point.no.tree, crs = 32188)
only.park = st_transform(x = only.park, crs = 32188)

# Get random points
set.seed(456)
rdm.n = 5
rdm.pt = st_sample(x = only.park, 
                   size = rdm.n)
mapView(only.park, 
        col.regions = c("green"))+
  mapView(rdm.pt)  # Random points


## ---- echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Add a buffer around the point we want to look at 
n = 10*2
min.buff = units::set_units(x = 0, value = "m")
max.buff = units::set_units(x = 100, value = "m")
buffer = st_buffer(x = c(selected.point,selected.point.no.tree), dist = units::set_units(x = 100, value = "m"))


mapView(only.park, 
        col.regions = c("green"))+
  mapView(c(selected.point,selected.point.no.tree), # Get the poitn that was added to be looked at 
          col.regions = c("red")) +
  mapView(buffer, col.regions = c("red")) 


## ---- echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Get random distance 
# rdm.disttmp = rexp(n = n, rate = 1)*20
# hist(rdm.disttmp)
set.seed(6543)
n = 10*2
min.buff = units::set_units(x = 0, value = "m")
max.buff = units::set_units(x = 100, value = "m")
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
# rmd.ptdf.sf = st_transform(rmd.ptdf.sf, crs = 32188)

set.seed(456)
mapView(only.park, col.regions = c("green"))+
  mapView(st_sample(only.park,5)) + 
  mapView(c(selected.point,selected.point.no.tree),col.regions = c("red")) + 
  mapView(buffer, col.regions = c("red")) + 
  mapview(rmd.ptdf.sf,col.regions = c("pink"))


## ---- echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(456)
# Add random points that are occupying the space of the polygon (grid )
rdm.pt = st_sample(x = only.park, 
                   size = 100,type ="hexagonal")

mapView(only.park, col.regions = c("green"))+
  mapView(rdm.pt) 


## ----simulate_sampling_function-------------------------------------------------------------------------------------------------------------------------------------------------------
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


## ----simulate_sampling_plots----------------------------------------------------------------------------------------------------------------------------------------------------------
par(mfrow=c(2,2), lwd = .3)
sample.mean.pop.est(x = x, n = n, sample.size = 1, ylim = 300)
sample.mean.pop.est(x = x, n = n, sample.size = 10, ylim = 300)
sample.mean.pop.est(x = x, n = n, sample.size = 50, ylim = 300)
sample.mean.pop.est(x = x, n = n, sample.size = 500, ylim = 300)


## ---- echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------
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


## ---- echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------
.printGen(rdm, target, i)

