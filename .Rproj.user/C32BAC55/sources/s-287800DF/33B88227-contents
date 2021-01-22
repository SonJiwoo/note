# 라이브러리 가져오기
library(tidyverse)
library(grid)
library(gridExtra)
theme_set(theme_minimal())

# Data
y <- c(93, 112, 122, 135, 122, 150, 118, 90, 124, 114)
# Sufficient Statistics
n <- length(y)
s2 <- var(y)
my <- mean(y)

# helper functions to sample from and evaluate
# scaled inverse chi-squared distribution
rsinvchisq <- function(n, nu, s2, ...) nu*s2 / rchisq(n , nu, ...)
dsinvchisq <- function(x, nu, s2){
  exp(log(nu/2)*nu/2 - lgamma(nu/2) + log(s2)/2*nu - log(x)*(nu/2+1) - (nu*s2/2)/x)
}

# sample 1000 random numbers from p(sigma2|y)
ns <- 1000
sigma2  <- rsinvchisq(ns, n-1, s2)
# sample 1000 from p(mu|sigma, y)
mu <- rnorm(n=length(sigma2), mean=my, sd=sqrt(sigma2/n))

# sample from predictive distribution p(ynew|y) for each draw of mu, sigma
sigma <- sqrt(sigma2)
ynew <- rnorm(ns, mu, sigma)

# grid range for mu, sigma and ynew
t1l <- c(90, 150)
t2l <- c(10, 60)
nl <- c(50, 185)
t1 <- seq(t1l[1], t1l[2], length.out = ns)
t2 <- seq(t2l[1], t2l[2], length.out = ns)
xynew <- seq(nl[1], nl[2], length.out = ns)

# compute exact marginal density of mu and sigma
# multiplication by 1./sqrt(s2/n) is due to the transformation of variable z=(x-mean(y))/sqrt(s2/n)
# the multiplication by 2*t2 is due to the transformation of variable z=t2^2,
pm <- dt((t1-my) / sqrt(s2/n), n-1) / sqrt(s2/n)
ps <- dsinvchisq(t2^2, n-1, s2) * 2*t2

# estimate the marginal density using samples and ad hoc Gaussian kernel approximation
pmk <- density(mu, adjust = 2, n = ns, from = t1l[1], to = t1l[2])$y
psk <- density(sigma, n = ns, from = t2l[1], to = t2l[2])$y

# compute the exact predictive density
# multiplication by 1./sqrt(s2/n) is due to the transformation of variable
p_new <- dt((xynew-my) / sqrt(s2*(1+1/n)), n-1) / sqrt(s2*(1+1/n))
                                                       
# Evaluate the joint density in a grid.
# Note that the following is not normalized, but for plotting contours it does not matter.
dfj <- data.frame(t1 = rep(t1, each = length(t2)), t2 = rep(t2, length(t1)))
dfj$z <- dsinvchisq(dfj$t2^2, n-1, s2) * 2*dfj$t2 * dnorm(dfj$t1, my, dfj$t2/sqrt(n))
# breaks for plotting the contours
cl <- seq(1e-5, max(dfj$z), length.out = 6)

# -----------------------------------------------------------------
# Demo 3.1. Visualize the joint and marginal densities
# Visualise the joint density and marginal densities of the posterior of normal distribution with unknown mean and variance.

# marginal density of mu
dfm <- data.frame(t1, Exact = pm, Empirical = pmk) %>% gather(grp, p, -t1)
margmu <- ggplot(dfm) +
  geom_line(aes(t1, p, color = grp)) +
  coord_cartesian(xlim = t1l) +
  labs(title = 'Marginal of mu', x = '', y = '') +
  scale_y_continuous(breaks = NULL) +
  theme(legend.background = element_blank(),
        legend.position = c(0.75, 0.8),
        legend.title = element_blank())
# marginal density of sigma
dfs <- data.frame(t2, Exact = ps, Empirical = psk) %>% gather(grp, p, -t2)
margsig <- ggplot(dfs) +
  geom_line(aes(t2, p, color = grp)) +
  coord_cartesian(xlim = t2l) +
  coord_flip() +
  labs(title = 'Marginal of sigma', x = '', y = '') +
  scale_y_continuous(breaks = NULL) +
  theme(legend.background = element_blank(),
        legend.position = c(0.75, 0.8),
        legend.title = element_blank())

# joint density
joint1labs <- c('Samples','Exact contour')
joint1 <- ggplot() +
  geom_point(data = data.frame(mu,sigma), aes(mu, sigma, col = '1'), size = 0.1) +
  geom_contour(data = dfj, aes(t1, t2, z = z, col = '2'), breaks = cl) +
  coord_cartesian(xlim = t1l,ylim = t2l) +
  labs(title = 'Joint posterior', x = '', y = '') +
  scale_y_continuous(labels = NULL) +
  scale_x_continuous(labels = NULL) +
  scale_color_manual(values=c('blue', 'black'), labels = joint1labs) +
  guides(color = guide_legend(nrow  = 1, override.aes = list(
    shape = c(16, NA), linetype = c(0, 1), size = c(2, 1)))) +
  theme(legend.background = element_blank(),
        legend.position = c(0.5, 0.9),
        legend.title = element_blank())

# combine the plots
# blank plot for combining the plots
bp <- grid.rect(gp = gpar(col = 'white'))
grid.arrange(joint1, margsig, margmu, bp, nrow = 2)


# -----------------------------------------------------------------
# Demo 3.2 Visualize Factored Distribution
# data frame for the conditional of mu and marginal of sigma
dfc <- data.frame(mu = t1, marg = rep(sigma[1], length(t1)),
                  cond = sigma[1] + dnorm(t1 ,my, sqrt(sigma2[1]/n)) * 100) %>%
  gather(grp, p, marg, cond)
# legend labels for the following plot
joint2labs <- c('Exact contour plot', 'Sample from joint post.',
                'Cond. distribution of mu', 'Sample from the marg. of sigma')
joint2 <- ggplot() +
  geom_contour(data = dfj, aes(t1, t2, z = z, col = '1'), breaks = cl) +
  geom_point(data = data.frame(m = mu[1], s = sigma[1]), aes(m , s, color = '2')) +
  geom_line(data = dfc, aes(mu, p, color = grp), linetype = 'dashed') +
  coord_cartesian(xlim = t1l,ylim = t2l) +
  labs(title = 'Joint posterior', x = '', y = '') +
  scale_x_continuous(labels = NULL) +
  scale_color_manual(values=c('blue', 'darkgreen','darkgreen','black'), labels = joint2labs) +
  guides(color = guide_legend(nrow  = 2, override.aes = list(
    shape = c(NA, 16, NA, NA), linetype = c(1, 0, 1, 1)))) +
  theme(legend.background = element_blank(),
        legend.position = c(0.5, 0.85),
        legend.title = element_blank())

margsig2 <- ggplot(data = data.frame(t2, ps)) +
  geom_line(aes(t2, ps), color = 'blue') +
  coord_cartesian(xlim = t2l) +
  coord_flip() +
  labs(title = 'Marginal of sigma', x = '', y = '') +
  scale_y_continuous(labels = NULL)

grid.arrange(joint2, margsig2, ncol = 2)


# -----------------------------------------------------------------
# Demo 3.3 Visualize the Marginal Distribution of Mu
# Visualise the marginal distribution of mu as a mixture of normals.
condpdfs <- sapply(t1, function(x) dnorm(x, my, sqrt(sigma2/n)))

# data frame of 25 first samples
dfm25 <- data.frame(t1, t(condpdfs[1:25,])) %>% gather(grp, p, -t1)
condmu <- ggplot(data = dfm25) +
  geom_line(aes(t1, p, group = grp), linetype = 'dashed') +
  labs(title = 'Conditional distribution of mu for first 25 samples', y = '', x = '') +
  scale_y_continuous(breaks = NULL)

# creata a plot of their mean
dfsam <- data.frame(t1, colMeans(condpdfs), pm) %>% gather(grp,p,-t1)
# labels
mulabs <- c('avg of sampled conds', 'exact marginal of mu')
meanmu <- ggplot(data = dfsam) +
  geom_line(aes(t1, p, size = grp, color = grp)) +
  labs(y = '', x = 'mu') +
  scale_y_continuous(breaks = NULL) +
  scale_size_manual(values = c(2, 0.8), labels = mulabs) +
  scale_color_manual(values = c('orange', 'black'), labels = mulabs) +
  theme(legend.position = c(0.8, 0.8),
        legend.background = element_blank(),
        legend.title = element_blank())
# combine the plots 
grid.arrange(condmu, meanmu, ncol = 1)

# -----------------------------------------------------------------
# Demo 3.4 Visualise Posterior Predictive Distirubtion
