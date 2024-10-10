rm(list=ls())
setwd('/Users/liuzhishuai/Desktop/DTR-MRL/code')
library('rgenoud')
library('survival')
library('speff2trial')
library('evd')
source('functions.R')

data(ACTG175)

Time <- c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800)
random_allocation <- rep(NA, length(Time))
MRL_IPW <- rep(NA, length(Time))
MRL_SIPW <- rep(NA, length(Time))
sd_IPW <- rep(NA, length(Time))
sd_SIPW <- rep(NA, length(Time))
#================Analysis 4=======================
data <- ACTG175

data$age <- data$age/max(data$age)
data$wtkg <- data$wtkg/max(data$wtkg)
data$karnof <- data$karnof/max(data$karnof)
data$cd40 <- data$cd40/max(data$cd40) 

attach(data)
A <- as.numeric(arms == 0 | arms == 3)
N <- nrow(data)
cens_rate <- sum(data$cens)/N
KM <- survfit(Surv(days, cens) ~ A)
plot(KM, lty = c(1,3))
legend(0.1, 0.3, c(1, 0), lty = c(1,3))
set.seed('1234567')
Y <- days + 0.1*rnorm(N, 0, 1)
prop <- sum(A == 1)/N
PS <- ifelse(A == 1, prop, 1 - prop)
delta <- cens

# Xclassify <- data.frame(karnof, drugs)
# Xclassify <- data.frame(wtkg, drugs, age) 
Xclassify <- data.frame(age, wtkg, karnof, cd40)

dim <- ncol(Xclassify) + 1

status <- 1 - cens
KM <- survfit(Surv(Y, status) ~ 1)
Sc <- rep(NA, N)
for (i in 1:N){
  temp <- which(Y[i] <= KM$time)[1]
  if (KM$surv[temp] == 0){
    temp <- temp - 1
  } 
  Sc[i] <- 1/KM$surv[temp]
}

for(i in 1:length(Time)){
  t <- Time[i]
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
  
  random_allocation[i] <- mean((Y - t)*(Y > t)*delta*Sc)/mean((Y > t)*delta*Sc)
  MRL_IPW[i] <- qValue(eta_opt)
  MRL_SIPW[i] <- qValue(Seta_opt)
  
  Variance1 <- ESTSTD(Y = Y, X.PS = matrix(1, N, 1), Xclassify = Xclassify, delta = delta, eta_hat = eta_opt,
                      treat = A, PS = PS, Sc = Sc, theta_hat = prop, t = t)
  SVariance1 <- Smooth_ESTSTD(Y = Y, X.PS = matrix(1, N, 1), Xclassify = Xclassify, delta = delta, eta_hat = Seta_opt,
                              treat = A, PS = PS, Sc = Sc, theta_hat = prop, t = t)
  
  sd_IPW[i] <- Variance1$sd
  sd_SIPW[i] <- SVariance1$sd
  
  Result <- rbind(random_allocation, MRL_IPW, sd_IPW, MRL_SIPW, sd_SIPW)
  
}
save(Result, file = '1-Result4.RData')

detach(data)