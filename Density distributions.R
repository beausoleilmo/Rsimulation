# Description  ------------------------------------------------------------
#### ### ### ## #### ### ### ## #### ### ### ## 
# Draw Density distributions 
# Created by Marc-Olivier Beausoleil
# 
# Why: 
  # To get a better sense of the distributions, draw the distributions based on their parameters 
# Requires:
# NOTES: 
# Reference : 
  # https://www.statology.org/r-guides/
#### ### ### ## #### ### ### ## #### ### ### ## 


# Binomial ----------------------------------------------------------------
# https://www.statology.org/plot-binomial-distribution-r/
# https://sphweb.bumc.bu.edu/otlt/mph-modules/bs/bs704_probability/bs704_probability7.html
success <- 0:25
plot(success, dbinom(x = success, size=20, prob=0.4),type='h', ylim = c(0,1)) 
plot(factor(success), dbinom(x = success, size=20, prob=0.5),type='h', ylim = c(0,1)) 
# Asks the question : what is the probability of getting k successes in a row given a certain probability of the object we are using to get that answer

# With size = 1 
# Example : fair coin (prob = 0.5) flipped once (size = 1)
plot(factor(0:1), dbinom(x = 0:1, size = 1, prob=0.5),type='h', ylim = c(0,1)) 
abline(h=seq(0,1, by = .1), lty = 3, lwd = .3)

# probability of getting  3 heads in a row (result = 1, three times in a row) with a fair coin (prob = 0.5), flipped 2 times (size = 3)
plot(factor(0:4), dbinom(x = 0:4, size = 3, prob=0.5),type='h', ylim = c(0,1)) 
abline(h=seq(0,1, by = .1), lty = 3, lwd = .3)
# Adding the calculation by hand : 
# The probability of getting 3 heads would be the same as having 1/2 chance of getting head 3 times in a row (1/2*1/2*1/2) = (1/2)^3
points(x = 4,(1/2)^3, pch = 19) 
# Here I added "4" although we never flip a sequence of 4 coins. Therefore, the probability is for sure going to be 0. 
# And from the plot, the probability of 4 tosses is actually 0. 

# Logistic distribution ---------------------------------------------------
# https://bookdown.org/roback/bookdown-BeyondMLR/ch-logreg.html 

curve(dlogis(x,location = 0, scale = 1), from=-10, to=10)
# hist(rlogis(n))
n = 1000
x = rnorm(n,0)
a=0
b=2
# c=-2
z = a + b*x #+ c*x^2
# location
m = 0
# scale 
s = 1

# this is the logit link 
pr = 1/(1+exp(-(z-m)/s))    
y = rbinom(n,1,pr)
col = "red"
plot(x,y)
# Why is a line (linear regression) is a problem for this? 
# Because, first the probabilities CAN'T be below 0 or above 1
lm.out = lm(y~x)
abline(lm.out, lty = 3)
# Logistic regression is a binomial regression with the "logistic" link function so in R (binomial(link = "logit"))
go=glm( y~x,#+I((x)^2),
        family="binomial")
newdata <- data.frame(x=seq(min(x), max(x),len=n))
newdata$y = predict(go,newdata, type="response") 
lines(x = newdata$x,
      y = newdata$y, 
      col = col,
      lwd = 2, 
      ylim = c(0,1))
######

set.seed(1234) 
n <- 500
b0 <- 0.5
b1 <- 1
x <- rnorm(n, mean = 10, sd = 2) 
ystar <- b0 + b1 * x + rlogis(n, location = -10)
y <- 1 * (ystar > 0)
mydat <- data.frame(x, y)
plot(mydat)

# Chi-square --------------------------------------------------------------
# https://www.statology.org/plot-chi-square-distribution-in-r/
curve(dchisq(x, df = 10), from = 0, to = 40) 


# Exponential  ------------------------------------------------------------
# https://www.statology.org/plot-exponential-distribution-in-r/
#plot PDF curves 
curve(dexp(x, rate = .5), from=0, to=10, col='blue')
curve(dexp(x, rate = 1), from=0, to=10, col='red', add=TRUE)
curve(dexp(x, rate = 1.5), from=0, to=10, col='purple', add=TRUE)
#add legend
legend(7, .5, legend=c("rate=.5", "rate=1", "rate=1.5"),
       col=c("blue", "red", "purple"), lty=1, cex=1.2)


# F-dsitribution ----------------------------------------------------------
# https://www.geo.fu-berlin.de/en/v/soga/Basics-of-statistics/Continous-Random-Variables/F-Distribution/F-Distribution-in-R/index.html
curve(df(x, df1 = 10, df2 = 20), from = 0, to = 4, n = 5000, col= 'pink', lwd=2)


# Normal ------------------------------------------------------------------
# https://stackoverflow.com/questions/9046664/how-to-use-the-function-curve-in-r-to-graph-a-normal-curve
# https://r-coder.com/normal-distribution-r/

curve(expr = dnorm(x = x, mean=10,sd=1), from = 5, to = 15, col="blue",main = "Density normal", ylim = c(0,1)) 

x <- seq(-5, 5, 0.1)

plot(x, dnorm(x, 0, 1), type = "l", lwd = 2, col = "blue", ylab = "", xlab = "x", ylim = c(0,1))
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
text(x = ub, y = p+0.04,labels = paste0(p*100,"%"),adj = 0,pos = 2)

# add the polygon to the right 
lb <- qnorm(1-p) # Lower bound
ub <- max(x)   # Upper bound
x2 <- seq(lb, max(x), length = 100) # New Grid
y <- dnorm(x2, 0, 1) # Density
polygon(c(lb, x2, ub), c(0,y,0), col = rgb(0, 0, 1, alpha = 0.5))
text(x = lb, y = p+0.04,labels = paste0(p*100,"%"),adj = 0,pos = 4)


# Log normal  -------------------------------------------------------------
# https://www.statology.org/plot-log-normal-distribution-r/
curve(dlnorm(x, meanlog=0, sdlog=1), from=0, to=10)


# Poisson -----------------------------------------------------------------
# https://www.statology.org/plot-poisson-distribution-r/

# The following reference is an excellent introduction to what are poisson distributions 
# https://bookdown.org/roback/bookdown-BeyondMLR/ch-poissonreg.html

#define range of "successes"
success <- 0:20
#create plot of probability mass function
plot(success, dpois(success, lambda=5), type='h')

set.seed(1)
n=200
x = rnorm(n)
beta0 = 2
beta1 = 0.5

log.mu = beta0 + beta1 *x
y = rpois(n, exp(log.mu))
plot(x,y, pch = 19)
# abline(exp(.5),exp(.3))

# poisson(link = "log")
glm.out = glm(y~x,family = poisson)
predProbs <- predict(glm.out,
                     newdata = data.frame(x=seq(min(x), max(x), length.out=100)), 
                     type="response")
lines(x = (seq(min(x), max(x), length.out=100)),
      y = (predProbs), col=2, lwd=2)

# glm.out$fitted.values
# abline(glm.out)
# lines(x,glm.out$fitted.values)



# t-distribution ----------------------------------------------------------
# https://www.statology.org/plot-t-distribution-r/
curve(dt(x, df=10), from=-4, to=4, col=col, lwd=lwd, ylim = c(0,1), ylab = "Probability", main = "t-distribution df = 10") 

# Uniform -----------------------------------------------------------------
# https://www.statology.org/uniform-distribution-r/
curve(dunif(x, 10,20), from=5, to=25) 



# All density distributions -----------------------------------------------
pdf("workshopXX-En/images/distributions_all.pdf",pointsize = 12, width = 9,height = 7)
# png("workshopXX-En/images/distributions_all.png",pointsize = 12, width = 9,height = 7, units = "in",res = 300)
par(mfrow = c(3,3), bg = NA)
# par(mfrow = c(1,1), bg = NA)
col = "black"
lwd=2
# Binom
plot(success, dbinom(x = success, size=20, prob=0.4),type='h', col = col, lwd=lwd, ylim = c(0,1), ylab = "Density", main = "Binomial, n = 20, p = .4")
curve(expr = dbinom(x = x, size=20, prob=0.6),col = "red", lwd=lwd, ylim = c(0,1), type = 'h',add = T)
legend("topright",legend = c("n = 20, p = .4","n = 20, p = .6"),lty = c(1,1), col =c("black","red"), lwd = 2)

# Logistic
curve(dlogis(x,location = 0, scale = 1), from=-10, to=10, col = col, lwd=lwd, ylim = c(0,1), ylab = "Density", main = "Logistic, l = 0, s = 1")

# Chi-sq
curve(dchisq(x, df = 10), from = 0, to = 40, col = col, lwd=lwd, ylim = c(0,1), ylab = "Density", main = "Chi-square, df = 10")
curve(dchisq(x, df = 4), from = 0, to = 40, col = "red", lwd=lwd, ylim = c(0,1), add = T)
legend("topright",legend = c("df = 10","df = 4"),lty = c(1,1), col =c("black","red"), lwd = 2)

# Exponential 
curve(dexp(x, rate = .5), from=0, to=10, col=col, lwd=lwd, ylim = c(0,1), ylab = "Density", main = "Exponential, rate = 0.5")
curve(dexp(x, rate = .2), from=0, to=10, col="red", lwd=lwd, ylim = c(0,1),add = T)
curve(dexp(x, rate = .8), from=0, to=10, col="blue", lwd=lwd, ylim = c(0,1),add = T)
legend("topright",legend = c("rate = 0.8","rate = 0.5","rate = 0.2"),lty = c(1,1,1), col =c("blue","black","red"), lwd = 2)

curve(df(x, df1 = 10, df2 = 20), from = 0, to = 4, n = 5000, col= col, lwd=lwd, ylim = c(0,1), ylab = "Density", main = "F-distribution, df1 = 10, df2 = 20")

# normal 
curve(expr = dnorm(x = x, mean=0,sd=1), from = -5, to = 5, col=col, lwd=lwd, ylim = c(0,1), ylab = "Density", main = "Normal, m = 0, sd = 1")
curve(expr = dnorm(x = x, mean=0,sd=2), from = -5, to = 5, col="red", lwd=lwd, ylim = c(0,1), add = T)
curve(expr = dnorm(x = x, mean=2,sd=1), from = -5, to = 5, col="green", lwd=lwd, ylim = c(0,1), add = T)
curve(dt(x, df=1), from=-5, to=5, col="blue", lty = 1, lwd=1, ylim = c(0,1), ylab = "Density", main = "t-distribution df = 10", add=TRUE)
legend("topright",legend = c("Normal, m = 0, sd = 1","Normal, m = 0, sd = 2","Normal, m = 2, sd = 1","t-distribution, df =1"),lty = c(1,1,1,1), col =c("black","red","green","blue"), lwd = 2)

curve(dlnorm(x, meanlog=0, sdlog=1), from=0, to=10, col=col, lwd=lwd, ylim = c(0,1), ylab = "Density", main = "Log-normal, m = 0 sd = 1")

plot(success, dpois(success, lambda=5), type='h', col=col, lwd=lwd, ylim = c(0,1), ylab = "Density", main = "Poisson lambda = 5")
# unifrom
curve(dunif(x, min = 10,max = 15), from=5, to=25, col=col, lwd=lwd, ylim = c(0,1), ylab = "Density", main = "Uniform min = 10, max = 15") 
curve(dunif(x, min = 16,max = 18), from=5, to=25, col="red", lwd=1, ylim = c(0,1), add = T) 
legend("topright",legend = c("min = 10, max = 15","min = 16, max = 18"),lty = c(1,1), col =c("black","red"), lwd = 2)
dev.off()





