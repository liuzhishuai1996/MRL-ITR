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
t <- 0
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
Y <- (days + 0.1*rnorm(N, 0, 1))
prop <- sum(A == 1)/N
PS <- ifelse(A == 1, prop, 1 - prop)
delta <- cens
proportion <- sum(Y > t)/N



# Xclassify <- data.frame(karnof, drugs)
# Xclassify <- data.frame(wtkg, drugs, age) 
Xclassify <- data.frame(age, wtkg, karnof, cd40)

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

########fit MRL########
data <- data.frame(Y = Y, X = rep(1, N), Z1 = karnof, Z2 = drugs, AZ1 = A*karnof, AZ2 = A*drugs, delta = delta, Sc = Sc)
para <- MRL(data = data, TIME = t, N2 = N, D1 = 0.01, D2 = 12, type = "proportional")

index <- which((Y[order(Y)] > t) == TRUE)[1]
TT <- Y[order(Y)]; TT <- c(0, TT) 



#AIPWE of Value-function(\tilde{V}(g,t))
qValue <- function(eta){
  treat_eta <- as.numeric(as.matrix(cbind(rep(1, N), Xclassify)) %*% eta > 0) 
  w <- (as.numeric(A)*treat_eta + (1 - as.numeric(A))*(1 - treat_eta))/PS
  MRL_hat <- exp(cbind(data$X, data$Z1, treat_eta*data$Z1, treat_eta*data$Z2) %*% para$est)
  m0 <- exp(cbind(data$X, data$Z1, treat_eta*data$Z1, treat_eta*data$Z2) %*% c(para$alpha[1], para$est[2:4]))
  mt <- exp(cbind(data$X, data$Z1, treat_eta*data$Z1, treat_eta*data$Z2) %*% c(para$alpha[index], para$est[2:4]))
  temp <- 0
  mu <- exp(data$X %*% t(para$alpha[1:index]) + cbind(data$Z1, treat_eta*data$Z1, treat_eta*data$Z2) %*% rbind(rep(para$est[2], index), rep(para$est[3], index), rep(para$est[4], index)))
  dT <- diff(TT[1:(index + 1)])
  temp <- (1/mu) %*% dT
  # for (i in 1:index){
  #   mu <- exp(cbind(data$X, data$Z1, treat_eta*data$Z1, treat_eta*data$Z2) %*% c(para$alpha[i], para$est[2:4]))
  #   temp <- temp + 1/mu*(TT[i + 1] - TT[i])
  # }
  St <- exp(- temp)*m0/mt
  temp1 <- (Y > t)*w*delta*(Y - t)*Sc + (1 - w)*MRL_hat*St
  Temp1 <- (Y > t)*w*delta*Sc + (1 - w)*St
  V1 <- mean(temp1)/mean(Temp1)
  return(V1)
}

########Smooth version of AIPWE of Value-function########
SqValue <- function(eta){
  c0 <- 4^(1/3)
  etaX <- as.matrix(cbind(rep(1, N), Xclassify)) %*% eta
  h <- c0*N^(- 1/3)*sd(etaX)
  Streat_eta <- pnorm(etaX/h, mean = 0, sd = 1)
  Sw <- (as.numeric(A)*Streat_eta + (1 - as.numeric(A))*(1 - Streat_eta))/PS
  treat_eta <- as.numeric(as.matrix(cbind(rep(1, N), Xclassify)) %*% eta > 0)
  MRL_hat <- exp(cbind(data$X, data$Z1, treat_eta*data$Z1, treat_eta*data$Z2) %*% para$est)
  m0 <- exp(cbind(data$X, data$Z1, treat_eta*data$Z1, treat_eta*data$Z2) %*% c(para$alpha[1], para$est[2:4]))
  mt <- exp(cbind(data$X, data$Z1, treat_eta*data$Z1, treat_eta*data$Z2) %*% c(para$alpha[index], para$est[2:4]))
  temp <- 0
  mu <- exp(data$X %*% t(para$alpha[1:index]) + cbind(data$Z1, treat_eta*data$Z1, treat_eta*data$Z2) %*% rbind(rep(para$est[2], index), rep(para$est[3], index), rep(para$est[4], index)))
  dT <- diff(TT[1:(index + 1)])
  temp <- (1/mu) %*% dT
  
  # tempp <- 0
  # for (i in 1:index){
  #   mu <- exp(cbind(data$X, data$Z1, treat_eta*data$Z1, treat_eta*data$Z2) %*% c(para$alpha[i], para$est[2:4]))
  #   tempp <- tempp + 1/mu*(TT[i + 1] - TT[i])
  # }
  
  St <- exp(- temp)*m0/mt
  temp2 <- (Y > t)*Sw*delta*(Y - t)*Sc + (1 - Sw)*MRL_hat*St
  Temp2 <- (Y > t)*Sw*delta*Sc + (1 - Sw)*St
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


# eta_opt <- c(0.05527386, 0.73219462, -0.60381546, -0.19725658, 0.23944221)
# Seta_opt <- c(0.1485504, 0.678992, -0.660294, -0.2201277,	0.1801621)
eta_opt <- c(-0.6520285, 0.5006694,	-0.1580004,	0.5205847, -0.1679774)
Seta_opt <- c(-0.25777162, 0.05251564, -0.74364366,	0.59244456,	-0.16370538)
# eta_opt <- c(0.6123424,	0.1855659, -0.5996196, -0.4597779, 0.1402236)
# Seta_opt <- c(0.70569786,	0.14504887,	-0.33060645, -0.60385831,	-0.08370109)
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

########compute variance and cover rate########
data <- data.frame(Y = Y, X = rep(1, N), Z1 = karnof, Z2 = drugs, AZ1 = A*karnof, AZ2 = A*drugs, delta = delta, Sc = Sc)
TIME_seq <- c(Y[order(Y)][Y[order(Y)] < t], t)
indi <- para$indi
alpha <- para$alpha
tempest <- para$est


Variance <- ESTSTD_AIPW(data = data, Y = Y, X.PS = matrix(1, N, 1), Xclassify = Xclassify, delta = delta, eta_hat = eta_opt,
                        treat = A, PS = PS, Sc = Sc, theta_hat = prop, t = t, alpha = alpha, tempest = tempest,
                        TIME_seq = TIME_seq, type = 'proportional')
SVariance <- Smooth_ESTSTD_AIPW(data = data, Y = Y, X.PS = matrix(1, N, 1), Xclassify = Xclassify, delta = delta, eta_hat = Seta_opt,
                         treat = A, PS = PS, Sc = Sc, theta_hat = prop, t = t, alpha = alpha, tempest = tempest,
                         TIME_seq = TIME_seq, type = 'proportional')

Variance$sd
SVariance$sd

detach(data)