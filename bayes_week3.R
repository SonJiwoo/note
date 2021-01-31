library(tidyverse)
library(gridExtra)
library(MASS)
library(reshape2)


### 1. Draw yourself Figure 7.1

# 초기 설정
inv <- solve
MU = matrix(c(50,50), ncol=1)
SIGMA = matrix(c(64,0,0,144), ncol=2)

# MVN pdf
calc.dmvn = Vectorize(function(a,b, mu=MU, sigma=SIGMA){
  y <- c(a,b)
  log.p <- (-nrow(mu)/2)*log(2*pi) - 0.5*log(det(sigma)) - 0.5*(t(y-mu) %*% inv(sigma) %*% (y-mu))
  exp(log.p)
})


# do it at once
allInOne <- function(corr){
  SIGMA = matrix(c(64,0,0,144), ncol=2)
  s1 <- sqrt(SIGMA[1,1]); s2 <- sqrt(SIGMA[2,2])
  SIGMA[1,2] <-  s1*s2*corr; SIGMA[2,1] <-  s1*s2*corr
  
  # MVN density function
  calc.dmvn = Vectorize(function(a,b, mu=MU, sigma=SIGMA){
    y <- c(a,b)
    log.p <- (-nrow(mu)/2)*log(2*pi) - 0.5*log(det(sigma)) - 0.5*(t(y-mu) %*% inv(sigma) %*% (y-mu))
    exp(log.p)
  })
  
  # sample
  sample = mvrnorm(n=30, mu=MU, Sigma=SIGMA)
  sample = data.frame(sample)
  colnames(sample) = c('y1','y2')
  # calculate density
  xLim = seq(20, 80, length=101)
  yLim = seq(20, 80, length=101)
  density.mvn <- outer(xLim, yLim, FUN=calc.dmvn)
  rownames(density.mvn) <- xLim
  colnames(density.mvn) <- yLim
  density.mvn <- melt(density.mvn)
  # graph
  density.mvn %>% 
    ggplot(aes(x=Var1, y=Var2)) +
    geom_tile(aes(fill=value, alpha=value)) +
    geom_contour(aes(z=value), color='white', size=0.1) +
    scale_fill_gradient(low='grey', high='black', guide=FALSE) +
    scale_alpha(guide=FALSE) +
    theme(legend.position='None') + theme_bw() +
    ggtitle(paste0('corr=',corr)) + xlab('y1') + ylab('y2') +
    geom_label(data=sample, aes(x=y1,y=y2,label='X'), size=1)
}
p1 <- allInOne(-0.5)
p2 <- allInOne(0)
p3 <- allInOne(0.5)
grid.arrange(p1,p2,p3, nrow=1) ((((((((()))))))))

### 2. Draw yourself Figure 7.2