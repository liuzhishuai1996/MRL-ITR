rm(list=ls())
setwd()
library('rgenoud')
library('survival')
library('speff2trial')
library('evd')
#======functions=======
source('functions.R')

data(ACTG175)
data <- ACTG175
attach(data)
A <- as.numeric(arms == 0)
N <- nrow(data)
cens_rate <- sum(data$cens)/N
KM <- survfit(Surv(days, cens) ~ A)
plot(KM, lty = c(1,3))
legend(0.1, 0.3, c(1, 0), lty = c(1,3)) 
Y <- days
PS <- rep(0.5, N)
delta <- cens
Xclassify <- symptom
dim <- ncol(Xclassify) + 1

status <- 1 - cens
KM <- survfit(Surv(days, status) ~ 1)
Sc <- rep(NA, N)
for (i in 1:N){
  temp <- which(days[i] <= KM$time)[1]
  if (KM$surv[temp] == 0){
    temp <- temp - 1
  } 
  Sc[i] <- 1/KM$surv[temp]
}
Sc[status] <- 0

Xclassify <- cbind(age, wtkg, karnof, cd40, cd80, homo, hemo, drugs, race, gender, str2, symptom)
NAME <- c('age', 'wtkg', 'karnof', 'cd40', 'cd80', 'homo', 'hemo', 'drugs', 'race', 'gender', 'str2', 'symptom')
dim <- ncol(Xclassify) + 1
DATA <- data.frame(Y, Xclassify, A)
p <- rep(NA, 12)
for (i in 1:12){
  lmfit <- lm(log(Y) ~ Xclassify[, i]*A, weights = Sc*cens, data = DATA)
  # lmfit <- lm(Y ~ Xclassify[, i]*A, weights = Sc*cens, data = DATA)
  results <- summary(lmfit)
  p[i] <- results$coefficients[4,4]
}
NAME[order(p)[1:5]]
p[order(p)[1:5]]