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
t <- 700
#================Analysis 4=======================
attach(ACTG175)
data <- ACTG175[(arms == 1 | arms == 2), ]
detach(ACTG175)
#scale the covariates
data$age <- data$age/max(data$age)
data$wtkg <- data$wtkg/max(data$wtkg)
data$karnof <- data$karnof/max(data$karnof)
data$cd40 <- data$cd40/max(data$cd40) 

attach(data)
A <- as.numeric(arms == 1)
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
proportion <- sum(Y > t)/N



Xclassify <- data.frame(karnof, drugs)
# Xclassify <- data.frame(wtkg, drugs, age) 
# Xclassify <- data.frame(age, wtkg, karnof, cd40)

dim <- ncol(Xclassify) + 1


#fit survival function for censoring
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


# eta_opt <- genoud(qValue, nvars = dim, max = TRUE, pop.size = 1000, default.domains = 1, print.level = 0, 
# optim.method = "Nelder-Mead", starting.values = rep(0, dim))$par
# eta_opt <- eta_opt/sqrt(sum(eta_opt^2))

# Seta_opt <- genoud(SqValue, nvars = dim, max = TRUE, pop.size = 1000, default.domains = 1, print.level = 0, 
# starting.values = rep(0, dim))$par

# Seta_opt <- Seta_opt/sqrt(sum(Seta_opt^2))

# eta_opt
# Seta_opt

# qValue(eta_opt)
# SqValue(Seta_opt)

mean((Y - t)*(Y > t)*delta*Sc)/mean((Y > t)*delta*Sc)

# eta_opt <- c(0.341, 0.636, -0.540, -0.404, -0.156)
# Seta_opt <- c(0.149, 0.679, -0.660, -0.220, 0.180)
# eta_opt <- c(0.584, -0.734, 0.345)
# Seta_opt <- c(0.584, -0.775, 0.240)
eta_opt <- c(0.317, -0.340, 0.885)
Seta_opt <- c(0.686, -0.724, 0.074)


# eta_opt <- c(-0.642, 0.496, -0.044, 0.494, -0.310)
# Seta_opt <- c(-0.501, 0.560, -0.335, 0.460, -0.334)

# eta_opt <- c(0.129, 0.198, -0.694, 0.245, -0.635)
# Seta_opt <- c(0.041, 0.213, -0.754, 0.328, -0.526)

treat_eta <- as.numeric(as.matrix(cbind(rep(1, N), Xclassify)) %*% eta_opt > 0)
treat_Seta <- as.numeric(as.matrix(cbind(rep(1, N), Xclassify)) %*% Seta_opt > 0)
sum(treat_eta == treat_Seta)/N

MAT <- matrix(0, 2, 2)
for (i in 1:N){
  if (treat_eta[i] == 0 & treat_Seta[i] == 0){
    MAT[1,1] <- MAT[1,1] + 1
  }else if(treat_eta[i] == 0 & treat_Seta[i] == 1){
    MAT[1,2] <- MAT[1,2] + 1
  }else if(treat_eta[i] == 1 & treat_Seta[i] == 0){
    MAT[2,1] <- MAT[2,1] + 1
  }else{
    MAT[2,2] <- MAT[2,2] + 1
  }
}
MAT

# Variance1 <- ESTSTD(Y = Y, X.PS = matrix(1, N, 1), Xclassify = Xclassify, delta = delta, eta_hat = eta_opt,
# treat = A, PS = PS, Sc = Sc, theta_hat = prop, t = t)
# SVariance1 <- Smooth_ESTSTD(Y = Y, X.PS = matrix(1, N, 1), Xclassify = Xclassify, delta = delta, eta_hat = Seta_opt,
# treat = A, PS = PS, Sc = Sc, theta_hat = prop, t = t)


Variance1$sd
SVariance1$sd

detach(data)


