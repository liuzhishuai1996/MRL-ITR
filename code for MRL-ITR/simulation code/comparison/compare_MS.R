rm(list=ls())
setwd()
library('rgenoud')
library('survival')
library('evd')
source('functions.R')
library(compiler)
enableJIT(3)

N <- 500
t <- 0
t0 <- 2
p <- 1; q <- - 3; c <- 0.5; d <- 1
C0 <- 0 
C1 <- 8.5 
D1 <- - 0.1
D2 <- 1
iteration <- 1000
ETA_IPW <- matrix(NA, iteration, 3)
SETA_IPW <- matrix(NA, iteration, 3)
ETA_AIPW <- matrix(NA, iteration, 3)
SETA_AIPW <- matrix(NA, iteration, 3)
regenerate <- 5
REGENERATE <- FALSE
ETA <- matrix(NA, iteration, 3)
SETA <- matrix(NA, iteration, 3)
censor_rate <- rep(NA, iteration)
Mhat <- rep(NA, iteration)
Mhat_opt <- rep(NA, iteration)
Mhat_smooth <- rep(NA, iteration)
Mhat_opt_smooth <- rep(NA, iteration)
correctPSmodel <- TRUE
set.seed(123)
N_test <- 100000
X1_test <- runif(N_test, min = -1, max = 1)
X2_test <- runif(N_test, min = -1, max = 1)
Xclassify_test <- data.frame(X1_test, X2_test)
if (t0 > 0.1){
  opt.treat <- as.numeric(as.matrix(cbind(rep(1, N_test), Xclassify_test)) %*% c(0, 1, - 1) > 0)
}else{
  opt.treat <- as.numeric(as.matrix(cbind(rep(1, N_test), Xclassify_test)) %*% c(0, - 1, 1) > 0)
}

time_test <- rep(NA, N_test)
U_test <- runif(N_test, min = 0, max = 1)
l_test <- opt.treat*(X1_test - X2_test)

time_test <- (log((c*l_test - p)*log(U_test) + exp(q + l_test*d)) - (q + l_test*d))/(p - c*l_test)
V_test_opt <- sum(time_test[time_test > t0] - t0)/sum(time_test > t0)

acc_IPW_test <- rep(NA, iteration)
acc_AIPW_test <- rep(NA, iteration)
Sacc_IPW_test <- rep(NA, iteration)
Sacc_AIPW_test <- rep(NA, iteration)
V_IPW_test <- rep(NA, iteration)
V_AIPW_test <- rep(NA, iteration)
SV_IPW_test <- rep(NA, iteration)
SV_AIPW_test <- rep(NA, iteration)



Iter <- 0
iter <- 1
while (iter <= iteration){
  REGENERATE <- FALSE
  Iter <- Iter + 1
  set.seed(1000 + Iter)
  
  cat('The', iter, 'th iteration', '\n')
  
  X1 <- runif(N, min = -1, max = 1)
  X2 <- runif(N, min = -1, max = 1)
  
  prob <- expit(X1 - 0.5*X2)
  treat <- rbinom(n = N, size = 1, prob = prob)
  
  T <- rep(NA, N)
  U <- runif(N, min = 0, max = 1)
  l <- treat*(X1 - X2)
  
  T <- (log((c*l - p)*log(U) + exp(q + l*d)) - (q + l*d))/(p - c*l)
  
  C <- runif(N, min = C0, max = C1)
  
  delta <- T < C 
  censor_rate[iter] <- 1 - sum(delta)/N 
  Y <- ifelse(T < C, T, C) 
  
  ########REJECTION########
  for (r in 1:regenerate){
    if (!delta[order(Y)][N - r + 1]){
      REGENERATE <- TRUE
      break
    }
  }
  if (REGENERATE){
    next
  }
  
  Xalpha <- data.frame(X1, X2)
  if (correctPSmodel){
    PS.mod <- glm(treat ~ . - 1, family = binomial(link = 'logit'), data = Xalpha) # correct specified
    X.PS <- Xalpha
  }else{
    PS.mod <- glm(treat ~ 1, family = binomial(link = 'logit')) #misspecified
    X.PS <- as.matrix(rep(1, length(treat)))
  }
  PS <- ifelse(treat == 1, PS.mod$fitted.values, 1 - PS.mod$fitted.values)
  theta_hat <- PS.mod$coefficients
  # PS <- ifelse(A == 1, prob, 1 - prob)
  
  status <- 1 - delta
  KM <- survfit(Surv(Y, status) ~ 1)
  Sc <- rep(NA, N)
  for (i in 1:N){
    temp <- which(Y[i] <= KM$time)[1]
    if (KM$surv[temp] == 0){
      temp <- temp - 1
    } 
    Sc[i] <- 1/KM$surv[temp]
  }
  
  Xclassify <- data.frame(X1, X2)
  dim <- ncol(Xclassify) + 1
  
  data <- data.frame(Y = Y, X = rep(1, N), Z1 = X1, Z2 = X2, AZ1 = treat*X1, AZ2 = treat*X2, delta = delta, C = C, T = T, Sc = Sc)
  para <- MRL(data = data, TIME = t, N2 = N, D1 = D1, D2 = D2, type = "proportional")
  
  
  index <- which((Y[order(Y)] > t) == TRUE)[1]
  TT <- Y[order(Y)]; TT <- c(0, TT)
  
  #IPWE of Value-function
  qValue.IPW <- function(eta){
    treat_eta <- as.numeric(as.matrix(cbind(rep(1, N), Xclassify)) %*% eta > 0) 
    w <- (as.numeric(treat)*treat_eta + (1 - as.numeric(treat))*(1 - treat_eta))/PS
    temp1 <- (Y > t)*w*delta*(Y - t)*Sc
    Temp1 <- (Y > t)*w*delta*Sc
    V1 <- mean(temp1)/mean(Temp1)
    return(V1)
  }
  
  #Smooth version of IPWE of Value-function
  SqValue.IPW <- function(eta){
    c0 <- 4^(1/3)
    etaX <- as.matrix(cbind(rep(1, N), Xclassify)) %*% eta
    h <- c0*N^(- 1/3)*sd(etaX)
    Streat_eta <- pnorm(etaX/h, mean = 0, sd = 1)
    Sw <- (as.numeric(treat)*Streat_eta + (1 - as.numeric(treat))*(1 - Streat_eta))/PS
    temp2 <- (Y > t)*Sw*delta*(Y - t)*Sc
    Temp2 <- (Y > t)*Sw*delta*Sc
    V2 <- mean(temp2)/mean(Temp2)
    return(V2)
  }
  
  
  
  ########AIPWE of Value-function########
  qValue.AIPW <- function(eta){
    treat_eta <- as.numeric(as.matrix(cbind(rep(1, N), Xclassify)) %*% eta > 0) 
    w <- (as.numeric(treat)*treat_eta + (1 - as.numeric(treat))*(1 - treat_eta))/PS
    MRL_hat <- exp(cbind(data$X, data$Z1, treat_eta*data$Z1, treat_eta*data$Z2) %*% para$est)
    m0 <- exp(cbind(data$X, data$Z1, treat_eta*data$Z1, treat_eta*data$Z2) %*% c(para$alpha[1], para$est[2:4]))
    mt <- exp(cbind(data$X, data$Z1, treat_eta*data$Z1, treat_eta*data$Z2) %*% c(para$alpha[index], para$est[2:4]))
    temp <- 0
    mu <- exp(data$X %*% t(para$alpha[1:index]) + cbind(data$Z1, treat_eta*data$Z1, treat_eta*data$Z2) %*% rbind(rep(para$est[2], index), rep(para$est[3], index), rep(para$est[4], index)))
    dT <- diff(TT[1:(index + 1)])
    temp <- (1/mu) %*% dT
    St <- exp(- temp)*m0/mt
    temp1 <- (Y > t)*w*delta*(Y - t)*Sc + (1 - w)*MRL_hat*St
    Temp1 <- (Y > t)*w*delta*Sc + (1 - w)*St
    V1 <- mean(temp1)/mean(Temp1)
    return(V1)
  }
  
  ########Smooth version of AIPWE of Value-function########
  SqValue.AIPW <- function(eta){
    c0 <- 4^(1/3)
    etaX <- as.matrix(cbind(rep(1, N), Xclassify)) %*% eta
    h <- c0*N^(- 1/3)*sd(etaX)
    Streat_eta <- pnorm(etaX/h, mean = 0, sd = 1)
    Sw <- (as.numeric(treat)*Streat_eta + (1 - as.numeric(treat))*(1 - Streat_eta))/PS
    treat_eta <- as.numeric(as.matrix(cbind(rep(1, N), Xclassify)) %*% eta > 0) 
    MRL_hat <- exp(cbind(data$X, data$Z1, treat_eta*data$Z1, treat_eta*data$Z2) %*% para$est)
    m0 <- exp(cbind(data$X, data$Z1, treat_eta*data$Z1, treat_eta*data$Z2) %*% c(para$alpha[1], para$est[2:4]))
    mt <- exp(cbind(data$X, data$Z1, treat_eta*data$Z1, treat_eta*data$Z2) %*% c(para$alpha[index], para$est[2:4]))
    temp <- 0
    mu <- exp(data$X %*% t(para$alpha[1:index]) + cbind(data$Z1, treat_eta*data$Z1, treat_eta*data$Z2) %*% rbind(rep(para$est[2], index), rep(para$est[3], index), rep(para$est[4], index)))
    dT <- diff(TT[1:(index + 1)])
    temp <- (1/mu) %*% dT
    
    St <- exp(- temp)*m0/mt
    temp2 <- (Y > t)*Sw*delta*(Y - t)*Sc + (1 - Sw)*MRL_hat*St
    Temp2 <- (Y > t)*Sw*delta*Sc + (1 - Sw)*St
    V2 <- mean(temp2)/mean(Temp2)
    return(V2)
  }
  
  eta_opt <- genoud(qValue.IPW, nvars = dim, max = TRUE, pop.size = 1000, default.domains = 1, print.level = 0, 
                    optim.method = "Nelder-Mead", starting.values = rep(0, dim))$par
  eta_opt.IPW <- eta_opt/sqrt(sum(eta_opt^2))
  
  Seta_opt <- genoud(SqValue.IPW, nvars = dim, max = TRUE, pop.size = 1000, default.domains = 1, print.level = 0, 
                     starting.values = rep(0, dim))$par
  
  Seta_opt.IPW <- Seta_opt/sqrt(sum(Seta_opt^2))
  
  
  eta_opt <- genoud(qValue.AIPW, nvars = dim, max = TRUE, pop.size = 1000, default.domains = 1, print.level = 0, 
                    optim.method = "Nelder-Mead", starting.values = rep(0, dim))$par
  eta_opt.AIPW <- eta_opt/sqrt(sum(eta_opt^2))
  
  Seta_opt <- genoud(SqValue.AIPW, nvars = dim, max = TRUE, pop.size = 1000, default.domains = 1, print.level = 0, 
                     starting.values = rep(0, dim))$par
  
  Seta_opt.AIPW <- Seta_opt/sqrt(sum(Seta_opt^2))
  
  
  
  
  
  ETA_IPW[iter, ] <- eta_opt.IPW
  ETA_AIPW[iter, ] <- Seta_opt.IPW
  SETA_IPW[iter, ] <- eta_opt.AIPW
  SETA_AIPW[iter, ] <- Seta_opt.AIPW
  
  treat_IPW_test <- as.numeric(as.matrix(cbind(rep(1, N_test), Xclassify_test)) %*% eta_opt.IPW > 0)
  treat_AIPW_test <- as.numeric(as.matrix(cbind(rep(1, N_test), Xclassify_test)) %*% Seta_opt.IPW > 0)
  Streat_IPW_test <- as.numeric(as.matrix(cbind(rep(1, N_test), Xclassify_test)) %*% eta_opt.AIPW > 0)
  Streat_AIPW_test <- as.numeric(as.matrix(cbind(rep(1, N_test), Xclassify_test)) %*% Seta_opt.AIPW > 0)
  
  l_IPW_test <- treat_IPW_test*(X1_test - X2_test)
  T_IPW_test <- (log((c*l_IPW_test - p)*log(U_test) + exp(q + l_IPW_test*d)) - (q + l_IPW_test*d))/(p - c*l_IPW_test)
  
  l_AIPW_test <- treat_AIPW_test*(X1_test - X2_test)
  T_AIPW_test <- (log((c*l_AIPW_test - p)*log(U_test) + exp(q + l_AIPW_test*d)) - (q + l_AIPW_test*d))/(p - c*l_AIPW_test)
  
  l_IPW_test <- Streat_IPW_test*(X1_test - X2_test)
  ST_IPW_test <- (log((c*l_IPW_test - p)*log(U_test) + exp(q + l_IPW_test*d)) - (q + l_IPW_test*d))/(p - c*l_IPW_test)
  
  l_AIPW_test <- Streat_AIPW_test*(X1_test - X2_test)
  ST_AIPW_test <- (log((c*l_AIPW_test - p)*log(U_test) + exp(q + l_AIPW_test*d)) - (q + l_AIPW_test*d))/(p - c*l_AIPW_test)
  
  acc_IPW_test[iter] <- sum(opt.treat == treat_IPW_test)/N_test
  acc_AIPW_test[iter] <- sum(opt.treat == treat_AIPW_test)/N_test
  Sacc_IPW_test[iter] <- sum(opt.treat == Streat_IPW_test)/N_test
  Sacc_AIPW_test[iter] <- sum(opt.treat == Streat_AIPW_test)/N_test
  
  V_IPW_test[iter] <- sum(T_IPW_test[T_IPW_test > t0] - t0)/sum(T_IPW_test > t0)
  V_AIPW_test[iter] <- sum(T_AIPW_test[T_AIPW_test > t0] - t0)/sum(T_AIPW_test > t0)
  SV_IPW_test[iter] <- sum(ST_IPW_test[ST_IPW_test > t0] - t0)/sum(ST_IPW_test > t0)
  SV_AIPW_test[iter] <- sum(ST_AIPW_test[ST_AIPW_test > t0] - t0)/sum(ST_AIPW_test > t0)
  iter <- iter + 1
}

#IPW
eta0_hat.IPW <- paste(round(mean(ETA_IPW[, 1]), 3), paste('(', round(sd(ETA_IPW[, 1]), 3), ')', sep = ''))
eta1_hat.IPW <- paste(round(mean(ETA_IPW[, 2]), 3), paste('(', round(sd(ETA_IPW[, 2]), 3), ')', sep = ''))
eta2_hat.IPW <- paste(round(mean(ETA_IPW[, 3]), 3), paste('(', round(sd(ETA_IPW[, 3]), 3), ')', sep = ''))
CA.IPW <- paste(round(mean(acc_IPW_test), 3), paste('(', round(sd(acc_IPW_test), 3), ')', sep = ''))

#AIPW
eta0_hat.AIPW <- paste(round(mean(ETA_AIPW[, 1]), 3), paste('(', round(sd(ETA_AIPW[, 1]), 3), ')', sep = ''))
eta1_hat.AIPW <- paste(round(mean(ETA_AIPW[, 2]), 3), paste('(', round(sd(ETA_AIPW[, 2]), 3), ')', sep = ''))
eta2_hat.AIPW <- paste(round(mean(ETA_AIPW[, 3]), 3), paste('(', round(sd(ETA_AIPW[, 3]), 3), ')', sep = ''))
CA.SIPW <- paste(round(mean(Sacc_IPW_test), 3), paste('(', round(sd(Sacc_IPW_test), 3), ')', sep = ''))
#SIPW
eta0_hat.SIPW <- paste(round(mean(SETA_IPW[, 1]), 3), paste('(', round(sd(SETA_IPW[, 1]), 3), ')', sep = ''))
eta1_hat.SIPW <- paste(round(mean(SETA_IPW[, 2]), 3), paste('(', round(sd(SETA_IPW[, 2]), 3), ')', sep = ''))
eta2_hat.SIPW <- paste(round(mean(SETA_IPW[, 3]), 3), paste('(', round(sd(SETA_IPW[, 3]), 3), ')', sep = ''))
CA.AIPW <- paste(round(mean(acc_AIPW_test), 3), paste('(', round(sd(acc_AIPW_test), 3), ')', sep = ''))
#SAIPW
eta0_hat.SAIPW <- paste(round(mean(SETA_AIPW[, 1]), 3), paste('(', round(sd(SETA_AIPW[, 1]), 3), ')', sep = ''))
eta1_hat.SAIPW <- paste(round(mean(SETA_AIPW[, 2]), 3), paste('(', round(sd(SETA_AIPW[, 2]), 3), ')', sep = ''))
eta2_hat.SAIPW <- paste(round(mean(SETA_AIPW[, 3]), 3), paste('(', round(sd(SETA_AIPW[, 3]), 3), ')', sep = ''))
CA.SAIPW <- paste(round(mean(Sacc_AIPW_test), 3), paste('(', round(sd(Sacc_AIPW_test), 3), ')', sep = ''))


MRL_IPW <- paste(round(mean(V_IPW_test), 3), paste('(', round(sd(V_IPW_test), 3), ')', sep = ''))
MRL_AIPW <- paste(round(mean(V_AIPW_test), 3), paste('(', round(sd(V_AIPW_test), 3), ')', sep = ''))
SMRL_IPW <- paste(round(mean(SV_IPW_test), 3), paste('(', round(sd(SV_IPW_test), 3), ')', sep = '')) 
SMRL_AIPW <- paste(round(mean(SV_AIPW_test), 3), paste('(', round(sd(SV_AIPW_test), 3), ')', sep = '')) 


CR <- round(mean(censor_rate), 3)
V_test_opt <- round(V_test_opt, 3)


if (correctPSmodel){
  PS <- 'T'
}else{
  PS <- 'F'
}

colNames <- c('PS', 'eta0_hat', 'eta1_hat', 'eta2_hat', 'MRL(eta_hat)', 'CA', 'MRL')
rowNames <- c('IPWE', 'SIPWE', 'AIPWE', 'SAIPWE')
row1 <- c(PS, eta0_hat.IPW, eta1_hat.IPW, eta2_hat.IPW, MRL_IPW, CA.IPW, V_test_opt)
row2 <- c(PS, eta0_hat.AIPW, eta1_hat.AIPW, eta2_hat.AIPW, MRL_AIPW, CA.AIPW, V_test_opt)
row3 <- c(PS, eta0_hat.SIPW, eta1_hat.SIPW, eta2_hat.SIPW, SMRL_IPW, CA.SIPW, V_test_opt)
row4 <- c(PS, eta0_hat.SAIPW, eta1_hat.SAIPW, eta2_hat.SAIPW, SMRL_AIPW, CA.SAIPW, V_test_opt)

results <- data.frame(rbind(row1, row2, row3, row4), row.names = rowNames)
colnames(results) <- colNames

fileName <- paste('time-varying regime', 'MS', N, PS, t, CR, sep = '-')
fileName <- paste(fileName, 'csv', sep = '.')
#fileName <- paste('/Users/liuzhishuai/Documents/硕士/研三上学期/DTR-MRL/Review/data_codes_survival/', fileName, sep = '')
write.table(results, file = fileName, sep = ',', row.names = TRUE, col.names = NA)
