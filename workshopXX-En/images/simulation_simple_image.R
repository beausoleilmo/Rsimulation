# Simulation image (cheat Sheet)
# par(bg = NA)
set.seed(1234)
n = 100
x = rnorm(n, mean = 0, sd = 1)
beta0 = 3; beta1 = .9
y = beta0 + beta1 * x + rnorm(n)
plot(x,y, pch = 19, cex = .2)
lm.out = lm(y~x)
abline(lm.out)

# Check the "Line" assumptions 
# L: Linearity
plot(x,y, pch = 19, cex = 1)
abline(lm.out)
# I: Independence
  # Assumed 
# N: Normality
hist(lm.out$residuals) 

# E: Equality of variance 
plot(lm.out$residuals~lm.out$fitted.values, cex = 1, pch =19)
abline(h=0)

summary(lm.out)
# plot(lm.out)
plot(x,y, pch = 19, cex = 1)
abline(lm.out)

# Reference 
# https://www.youtube.com/watch?v=zQMu3ynel2s&t=634s