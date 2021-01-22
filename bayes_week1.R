# A First Course in Bayesian Statistical Methods
library(tidyverse)

## Exercise 3.3
####(a)
round(qgamma(c(0.025, 0.975), 237, 20),2)
round(qgamma(c(0.025, 0.975), 125, 14),2)

####(b)
n0 = 1:50
data.frame(n0, exp=(12*n0+113)/(n0+13)) %>% 
  ggplot(aes(x=n0, y=exp)) + 
  geom_point() + ylab(NULL) +
  ggtitle(bquote("Posterior Expectation of" ~ theta[B])) +
  geom_hline(yintercept = 237/20, color='red') +
  annotate("text", x=20, y=237/20-0.2,
           label=bquote("Criterion:Posterior Expectation of" ~ theta[A] ~"↗"), size=4)
# theta_B의 사후기댓값이 theta_A의 사후기댓값(11.83)과 비슷해지기 위해서는
# theta_B에 대한 사전믿음 수준이 11.83에 가깝고 그 강도도 더 강해야 할 것이다.

####(c)
# 집단A와 집단B가 사전기댓값이 12로 같기는 하지만, evidence 차이가 워낙 커서
# 두 집단은 독립이라고 봐도 될 정도인 것 같다.

#--------------------------------------------------------------------------------#

# Bayesian Data Analysis

## Exercise 1.1
####(a)
y <- seq(-5, 5, 0.05)
data.frame(y) %>% 
  ggplot(aes(x=y, 0.5*dnorm(y, mean=1, sd=2) + 0.5*dnorm(y,mean=2, sd=2))) +
  geom_line(size=1) + ylab(NULL)

####(b)
theta_1 = dnorm(x=1, mean=1, sd=2) # when theta=1
theta_2 = dnorm(x=1, mean=2, sd=2) # when theta=2
ans = theta_1 / (theta_1 + theta_2)
print(ans)

####(c)
# 분산이 커지면 평균이 1인지 2인지가 중요해지지 않는다.
# 그러므로 theta가 1일 확률과 2일 사후확률은 0.5로 점점 같아진다.
# 분산이 작아지면 데이터(y)의 값에 사후확률이 극단적으로 좌지우지된다.

## Exercise 2-13.
df = data.frame(year=seq(1976,1985),
           fatal.accident=c(24,25,31,31,22,21,26,20,16,22),
           pass.death=c(734,516,754,877,814,362,764,809,223,1066),
           death.rate=c(.19,.12,.15,.16,.14,.06,.13,.13,.03,.15))
print(df)

####(a)
accident.sum = sum(df[,"fatal.accident"])
year.cnt = length(df[,'year'])

sample.size = 1000
theta = rgamma(sample.size, accident.sum)/year.cnt
sample = rpois(sample.size, theta)
pred.interval = sort(sample)[c(25, 976)]
print(pred.interval) # 답: (14, 35) (물론 할 때마다 조금씩 다르게 나옴)

####(b)
df <- df %>% 
  mutate(pass.mile = round(pass.death/death.rate*1e8, -6))
print(df)

mile.sum = sum(df[,'pass.mile'])
theta = rgamma(sample.size, accident.sum)/mile.sum
sample = rpois(sample.size, theta*8e11) 
pred.interval = sort(sample)[c(25, 976)]
print(pred.interval) # 답: (21, 47) (물론 할 때마다 조금씩 다르게 나옴)

####(c)
death.sum = sum(df[,'pass.death'])
theta = rgamma(sample.size, death.sum)/year.cnt
sample = rpois(sample.size, theta)
pred.interval = sort(sample)[c(25, 976)]
print(pred.interval) # 답: (642, 745) (물론 할 때마다 조금씩 다르게 나옴)

####(d)
theta = rgamma(sample.size, death.sum)/mile.sum
sample = rpois(sample.size, theta*8e11)
pred.interval = sort(sample)[c(25, 976)]
print(pred.interval) # 답: (904, 1034) (물론 할 때마다 조금씩 다르게 나옴)
