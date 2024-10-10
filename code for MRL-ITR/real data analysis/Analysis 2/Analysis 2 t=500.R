rm(list=ls())
setwd()
library('rgenoud')
library('survival')
library('speff2trial')
library('evd')
#======functions=======
source('functions.R')

#data
data(ACTG175)

#pre-specified t
t <- 500
#================Analysis 4=======================
data <- ACTG175
#scaling
data$age <- data$age/max(data$age)
data$wtkg <- data$wtkg/max(data$wtkg)
data$karnof <- data$karnof/max(data$karnof)
data$cd40 <- data$cd40/max(data$cd40)

attach(data)
A <- as.numeric(arms == 0)
N <- nrow(data)
cens_rate <- sum(data$cens)/N
KM <- survfit(Surv(days, cens) ~ A)
plot(KM, lty = c(1,3))
legend(0.1, 0.3, c(1, 0), lty = c(1,3))
Y <- days + 0.1*rnorm(N, 0, 1)
prop <- sum(A == 1)/N
PS <- ifelse(A == 1, prop, 1 - prop)
delta <- cens
proportion <- sum(Y > t)/N


Xclassify <- data.frame(karnof)
# Xclassify <- data.frame(age, wtkg, karnof, cd40)
dim <- ncol(Xclassify) + 1


#fit survival function for censoring
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


#IPWE of Value-function(\tilde{V}(g,t))
qValue <- function(eta){
  treat_eta <- as.numeric(as.matrix(cbind(rep(1, N), Xclassify)) %*% eta > 0) 
  w <- (as.numeric(A)*treat_eta + (1 - as.numeric(A))*(1 - treat_eta))/PS
  temp1 <- (Y > t)*w*delta*(Y - t)*Sc
  Temp1 <- (Y > t)*w*delta*Sc
  V1 <- mean(temp1)/mean(Temp1)
  return(V1)
}

#Smooth version of IPWE of Value-function
SqValue <- function(eta){
  c0 <- 4^(1/3)
  etaX <- as.matrix(cbind(rep(1, N), Xclassify)) %*% eta
  h <- c0*N^(- 1/3)*sd(etaX)
  Streat_eta <- pnorm(etaX/h, mean = 0, sd = 1)
  Sw <- (as.numeric(A)*Streat_eta + (1 - as.numeric(A))*(1 - Streat_eta))/PS
  temp2 <- (Y > t)*Sw*delta*(Y - t)*Sc
  Temp2 <- (Y > t)*Sw*delta*Sc
  V2 <- mean(temp2)/mean(Temp2)
  return(V2)
}

eta_opt <- genoud(qValue, nvars = dim, max = TRUE, pop.size = 1000, default.domains = 1, print.level = 0, 
                  optim.method = "Nelder-Mead", starting.values = rep(0, dim))$par
eta_opt <- eta_opt/sqrt(sum(eta_opt^2))

Seta_opt <- genoud(SqValue, nvars = dim, max = TRUE, pop.size = 1000, default.domains = 1, print.level = 0, 
                   starting.values = rep(0, dim))$par

Seta_opt <- Seta_opt/sqrt(sum(Seta_opt^2))

eta_opt
Seta_opt

qValue(eta_opt)
SqValue(Seta_opt)

mean((Y - t)*(Y > t)*delta*Sc)/mean((Y > t)*delta*Sc)


#compute variance
Variance1 <- ESTSTD(Y = Y, X.PS = matrix(1, N, 1), Xclassify = Xclassify, delta = delta, eta_hat = eta_opt,
treat = A, PS = PS, Sc = Sc, theta_hat = prop, t = t)
SVariance1 <- Smooth_ESTSTD(Y = Y, X.PS = matrix(1, N, 1), Xclassify = Xclassify, delta = delta, eta_hat = Seta_opt,
treat = A, PS = PS, Sc = Sc, theta_hat = prop, t = t)

Variance1$sd
SVariance1$sd

detach(data)