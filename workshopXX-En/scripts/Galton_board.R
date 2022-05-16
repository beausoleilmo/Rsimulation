library(ggplot2)
library(ggExtra)
library(dplyr)

set.seed(1234567)
n = 500 # changes the sd! (number of pegs that got hit, will change the outcome in X)
rep.n = 1000 # nb.balls balls
sd = (sqrt(n)-1) # Get the standard deviation 

# Get the balls to hit a peg and move right (+1) or left (-1) 
pos.ball <- function(n) {
  rb = rbinom(n, size = 1, .5)
  pos.final = sum(ifelse(rb==1,1,-1))
  return(pos.final)  
}
# Custom density normal distribution
dn <- function(x,mean=0,sd=2) {
  dnorm(x,mean,sd)/max(dnorm(x,mean,sd))
}

# Make multiple balls bounce!
pop.ball = replicate(n = rep.n, expr = pos.ball(n = n), simplify = TRUE)

# Check if normal 
# hist(pop.ball, breaks = 40)
shapiro.test(pop.ball)
# qqnorm(pop.ball, main='Normal')
# qqline(pop.ball)
ks.test(x = pop.ball, y = 'pnorm')

df = data.frame(x = pop.ball)
yheight <- max(dplyr::count(df,x)["n"]) 

binwidth = 2.1
dotsize=1

bk = seq(0, 1, 1/yheight)
lb = seq(0, yheight, by = 1)
select.every = 5
select.odd = ifelse(1:length(bk) %% select.every == 1,TRUE,FALSE)

galtonqunincunx = ggplot() +
  geom_dotplot(data = df, aes(x),
               binwidth=binwidth, method="histodot", 
               dotsize = dotsize,
               fill = "white", stroke = 2) +
  coord_fixed(ratio=binwidth*dotsize*yheight)+
  theme_bw() +
  theme(panel.background=element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(1,0,1,0), "cm"),
        axis.ticks = element_line(colour = "black"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16),
        axis.line = element_line(colour = "black"))+ 
  # scale_x_continuous(breaks = seq(  range(pop.ball)[1], range(pop.ball)[2],5)) +
  scale_y_continuous(limits=c(0, 1), expand = c(0, 0), 
                     breaks = bk[select.odd], 
                     labels = lb[select.odd]) +
  removeGridX()+
  stat_function(fun = dn, args = list(mean = 0, sd = sd), geom = "line")+
  # labs(x=NULL, y=NULL) + 
  labs(title = "Galton quincunx (board)", x = "Sum of binomial outcome", y = "Count");galtonqunincunx
