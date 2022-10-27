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