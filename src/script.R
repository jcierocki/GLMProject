---
title: "GLM"
author: "JEH"
date: "2022/2/28"
---


## Viewing

setwd(" ")
dat <- read.table("GLM_dataset_16_victims.txt", header = TRUE)
print(dat)
summary(dat)


## Data entry

black <- c(119,16,12,7,3,2,0)
white <- c(1070,60,14,4,0,0,1)
resp <- c(rep(0:6,times=black), rep(0:6,times=white))
race <- factor(c(rep("black", sum(black)), rep("white", sum(white))), levels = c("white","black"))
victim <- data.frame(resp, race)


## exploring

table(race)
with(victim, tapply(resp, race, mean))
with(victim, tapply(resp, race, var))

table(resp, race)


## Poisson model

poisson_model <- glm(resp ~ race, data = victim, family = poisson)
summary(poisson_model)


## Ratio of sample means

exp(coef(poisson_model)[2])
mean(victim$resp[victim$race == "black"])/mean(victim$resp[victim$race == "white"])


# Fitted counts for Poisson GLM

fitted_means <- exp(predict(poisson_model, newdata = data.frame(race = c("white","black"))))
fitted_means

fitted_W <- dpois(0:6,lambda = fitted_means[1]) * sum(victim$race =="white") 
fitted_B <- dpois(0:6,lambda = fitted_means[2]) * sum(victim$race =="black") 
data.frame(Response = 0:6,BlackObs = black, BlackFit = round(fitted_B,1), WhiteObs = white, WhiteFit = round(fitted_W,1))


## rootogram

countreg::rootogram(poisson_model)


## negative binomial model

library(MASS)
neg_bin_model <- glm.nb(resp ~ race, data = victim)
summary(neg_bin_model)


# fitted counts for Negative Binomial GLM

fitted_means_nb <- exp(predict(poisson_model, newdata = data.frame(race = c("white","black"))))
fitted_means_nb 


## dispersion parameter

neg_bin_model$theta


## estimated variance for the count

fitted_means_nb + fitted_means_nb^2 * (1/neg_bin_model$theta)


## simulating the same number of observations 

library(magrittr)
op <- par(mfrow = c(1,2))
set.seed(1)

victim$resp %>% `[`(victim$race == "white") %>% 
  table() %>% barplot(main = "Observed White")
rnbinom(n = 1149, size = neg_bin_model$theta, mu = exp(coef(neg_bin_model)[1])) %>% 
  table() %>%  barplot(main = "Simulated White")

victim$resp %>% `[`(victim$race == "black") %>% 
  table() %>% barplot(main = "Observed Black")
rnbinom(n = 159, size = neg_bin_model$theta, mu = exp(sum(coef(neg_bin_model)))) %>% 
  table() %>%  barplot(main = "Simulated Black")
par(op)

