# Use base R here 
# https://stackoverflow.com/questions/71052975/how-to-plot-histograms-in-base-r-in-the-margin-of-a-plot?noredirect=1#comment125604177_71052975 
source(file = "workshopXX-En/scripts/barplot2.R")
set.seed(123)
# par(mfrow = c(1,2), mar =c(4,4,3,3), cex = 1.4)
n = 250
x.1 = rnorm(n, mean = 15, sd = 5)
y.1 = rnorm(n, mean = 5, sd = 2)

x.2 = rnorm(n, mean = 15, sd = 1)
y.2 = rnorm(n, mean = 5, sd = .5)

size.line = .8
col.line = alpha("black",.5)
df.1 <- data.frame(x = x.1, y = y.1)
df.2 <- data.frame(x = x.2, y = y.2)
data = df.1
layout(
  matrix(c(2,0,
           1,3), byrow=TRUE, nrow=2),
  widths = c(4,1), heights = c(1,4)
)
# main plot
par(mar = c(4, 5, 0, 0) + 0.1)
plot(y ~ x, data = data, pch = 16)
usr <- par("usr")
# top margin
par(mar = c(0, 5, 0, 0) + 0.1)
hx <- hist(data$x, plot = FALSE)
hy <- hist(data$y, plot = FALSE)
plot(NA, type = "n", xlim = usr[1:2], ylim = c(0, max(c(hx$counts, hy$counts))),
     xaxs = "i", yaxs = "i", xaxt = "n", ylab = "Counts", frame = FALSE)
barplot2(hx$mids, hx$counts, space = 0, horiz = FALSE, col = "gray")
# right margin
par(mar = c(4, 0, 0, 0) + 0.1)
plot(NA, type = "n", xlim = c(0, max(c(hx$counts, hy$counts))), ylim = usr[3:4],
     xaxs = "i", yaxs = "i", yaxt = "n", xlab = "Counts", frame = FALSE)
barplot2(hy$mids, hy$counts, space = 0, horiz = TRUE, col = "gray")


