rm(list=ls())
setwd()
library('rgenoud')
library('survival')
library('evd')
#======functins=======
source('functions.R')
#=======parameters=========
N <- 500
t <- 2
p <- 1; q <- - 3; c <- 0.5; d <- 1
C0 <- 0 
# C1 <- 17
C1 <- 8.5 
ETA <- matrix(NA, iteration, 3)
SETA <- matrix(NA, iteration, 3)
regenerate <- 1
REGENERATE <- FALSE
iteration <- 1000
ETA <- matrix(NA, iteration, 3)
SETA <- matrix(NA, iteration, 3)
censor_rate <- rep(NA, iteration)
Mhat <- rep(NA, iteration)
Mhat_opt <- rep(NA, iteration)
Std <- rep(NA, iteration)
cover <- rep(NA, iteration)
Mhat_smooth <- rep(NA, iteration)
Mhat_opt_smooth <- rep(NA, iteration)
Std_smooth <- rep(NA, iteration)
cover_smooth <- rep(NA, iteration)

Mhat1 <- rep(NA, iteration)
Mhat_opt1 <- rep(NA, iteration)
Std1 <- rep(NA, iteration)
cover1 <- rep(NA, iteration)
Mhat_smooth1 <- rep(NA, iteration)
Mhat_opt_smooth1 <- rep(NA, iteration)
Std_smooth1 <- rep(NA, iteration)
cover_smooth1 <- rep(NA, iteration)

Std2 <- rep(NA, iteration)
cover2 <- rep(NA, iteration)
Std_smooth2 <- rep(NA, iteration)
cover_smooth2 <- rep(NA, iteration)

correctPSmodel <- FALSE
#======generate testing set======
set.seed(123)
N_test <- 100000
X1_test <- runif(N_test, min = -1, max = 1)
X2_test <- runif(N_test, min = -1, max = 1)
Xclassify_test <- data.frame(X1_test, X2_test)
if (t > 0.1){
  opt.treat <- as.numeric(as.matrix(cbind(rep(1, N_test), Xclassify_test)) %*% c(0, 1, - 1) > 0)
}else{
  opt.treat <- as.numeric(as.matrix(cbind(rep(1, N_test), Xclassify_test)) %*% c(0, - 1, 1) > 0)
}

time_test <- rep(NA, N_test)
U_test <- runif(N_test, min = 0, max = 1)
l_test <- opt.treat*(X1_test - X2_test)

time_test <- (log((c*l_test - p)*log(U_test) + exp(q + l_test*d)) - (q + l_test*d))/(p - c*l_test)
V_test_opt <- sum(time_test[time_test > t] - t)/sum(time_test > t)

acc_test <- rep(NA, iteration)
Sacc_test <- rep(NA, iteration)
V_test <- rep(NA, iteration)
SV_test <- rep(NA, iteration)

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
  
  #IPWE of Value-function
  qValue <- function(eta){
    treat_eta <- as.numeric(as.matrix(cbind(rep(1, N), Xclassify)) %*% eta > 0)
    w <- (as.numeric(treat)*treat_eta + (1 - as.numeric(treat))*(1 - treat_eta))/PS
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
    Sw <- (as.numeric(treat)*Streat_eta + (1 - as.numeric(treat))*(1 - Streat_eta))/PS
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
  
  ETA[iter, ] <- eta_opt
  SETA[iter, ] <- Seta_opt
  
  ########Variance and cover rate########
  #\hat{M}_{IPWE}(\hat{eta_opt})
  Mhat[iter] <- qValue(eta_opt)
  Variance <- ESTSTD(Y = Y, X.PS = X.PS, Xclassify = Xclassify, delta = delta, eta_hat = c(0, 1, - 1),
                     treat = treat, PS = PS, Sc = Sc, theta_hat = theta_hat, t = t)
  Std[iter] <- Variance$sd
  Mhat_opt[iter] <- qValue(c(0 , 1, - 1))
  cover[iter] <- (V_test_opt < (Mhat_opt[iter] + 1.96*Variance$sd))&(V_test_opt > (Mhat_opt[iter] - 1.96*Variance$sd))
  
  Mhat_smooth[iter] <- SqValue(Seta_opt)
  SVariance <- Smooth_ESTSTD(Y = Y, X.PS = X.PS, Xclassify = Xclassify, delta = delta, eta_hat = c(0 , 1, - 1),
                             treat = treat, PS = PS, Sc = Sc, theta_hat = theta_hat, t = t)
  Std_smooth[iter] <- SVariance$sd
  Mhat_opt_smooth[iter] <- SqValue(c(0 , 1, - 1))
  cover_smooth[iter] <- (V_test_opt < (Mhat_opt_smooth[iter] + 1.96*SVariance$sd))&(V_test_opt > (Mhat_opt_smooth[iter] - 1.96*SVariance$sd))
  
  ########Variance and cover rate 1########
  Variance1 <- ESTSTD(Y = Y, X.PS = X.PS, Xclassify = Xclassify, delta = delta, eta_hat = eta_opt,
                      treat = treat, PS = PS, Sc = Sc, theta_hat = theta_hat, t = t)
  Std1[iter] <- Variance1$sd
  Mhat_opt1[iter] <- qValue(eta_opt)
  cover1[iter] <- (V_test_opt < (Mhat_opt1[iter] + 1.96*Variance1$sd))&(V_test_opt > (Mhat_opt1[iter] - 1.96*Variance1$sd))
  
  SVariance1 <- Smooth_ESTSTD(Y = Y, X.PS = X.PS, Xclassify = Xclassify, delta = delta, eta_hat = Seta_opt,
                              treat = treat, PS = PS, Sc = Sc, theta_hat = theta_hat, t = t)
  Std_smooth1[iter] <- SVariance1$sd
  Mhat_opt_smooth1[iter] <- SqValue(Seta_opt)
  cover_smooth1[iter] <- (V_test_opt < (Mhat_opt_smooth1[iter] + 1.96*SVariance1$sd))&(V_test_opt > (Mhat_opt_smooth1[iter] - 1.96*SVariance1$sd))
  
  ########bootstrap Variance and cover rate########
  Variance2 <- bootstrap_sd(Y = Y, Xclassify = Xclassify, delta = delta, treat = treat, Sc = Sc, PS = PS, w = w, eta_hat = eta_opt)
  Std2[iter] <- Variance2$boot_std1
  cover2[iter] <-  (V_test_opt < (Mhat_opt1[iter] + 1.96*Variance2$boot_std1))&(V_test_opt > (Mhat_opt1[iter] - 1.96*Variance2$boot_std1))
  
  Std_smooth2[iter] <- Variance2$boot_std2
  cover_smooth2[iter] <- (V_test_opt < (Mhat_opt_smooth1[iter] + 1.96*Variance2$boot_std2))&(V_test_opt > (Mhat_opt_smooth1[iter] - 1.96*Variance2$boot_std2))
  
  ########Performance########
  treat_test <- as.numeric(as.matrix(cbind(rep(1, N_test), Xclassify_test)) %*% eta_opt > 0)
  Streat_test <- as.numeric(as.matrix(cbind(rep(1, N_test), Xclassify_test)) %*% Seta_opt > 0)
  
  l_test <- treat_test*(X1_test - X2_test)
  T_test <- (log((c*l_test - p)*log(U_test) + exp(q + l_test*d)) - (q + l_test*d))/(p - c*l_test)
  
  l_test <- Streat_test*(X1_test - X2_test)
  ST_test <- (log((c*l_test - p)*log(U_test) + exp(q + l_test*d)) - (q + l_test*d))/(p - c*l_test)
  
  
  acc_test[iter] <- sum(opt.treat == treat_test)/N_test
  Sacc_test[iter] <- sum(opt.treat == Streat_test)/N_test
  V_test[iter] <- sum(T_test[T_test > t] - t)/sum(T_test > t)
  SV_test[iter] <- sum(ST_test[ST_test > t] - t)/sum(ST_test > t)
  iter <- iter + 1
}

########Results########
eta0_hat <- paste(round(mean(ETA[, 1]), 3), paste('(', round(sd(ETA[, 1]), 3), ')', sep = ''))
eta1_hat <- paste(round(mean(ETA[, 2]), 3), paste('(', round(sd(ETA[, 2]), 3), ')', sep = ''))
eta2_hat <- paste(round(mean(ETA[, 3]), 3), paste('(', round(sd(ETA[, 3]), 3), ')', sep = ''))

Seta0_hat <- paste(round(mean(SETA[, 1]), 3), paste('(', round(sd(SETA[, 1]), 3), ')', sep = ''))
Seta1_hat <- paste(round(mean(SETA[, 2]), 3), paste('(', round(sd(SETA[, 2]), 3), ')', sep = ''))
Seta2_hat <- paste(round(mean(SETA[, 3]), 3), paste('(', round(sd(SETA[, 3]), 3), ')', sep = ''))

CA <- paste(round(mean(acc_test), 3), paste('(', round(sd(acc_test), 3), ')', sep = ''))
SCA <- paste(round(mean(Sacc_test), 3), paste('(', round(sd(Sacc_test), 3), ')', sep = ''))

MRL <- paste(round(mean(V_test), 3), paste('(', round(sd(V_test), 3), ')', sep = ''))
SMRL <- paste(round(mean(SV_test), 3), paste('(', round(sd(SV_test), 3), ')', sep = '')) 

MRL_hat <- paste(round(mean(Mhat), 3), paste('(', round(sd(Mhat), 3), ')', sep = ''))
SMRL_hat <- paste(round(mean(Mhat_smooth), 3), paste('(', round(sd(Mhat_smooth), 3), ')', sep = ''))
if (correctPSmodel){
  PS <- 'T'
}else{
  PS <- 'F'
}
CP <- mean(cover)
SCP <- mean(cover_smooth)
SE <- round(mean(Std), 3)
SSE <- round(mean(Std_smooth), 3)

CP1 <- mean(cover1)
SCP1 <- mean(cover_smooth1)
SE1 <- round(mean(Std1), 3)
SSE1 <- round(mean(Std_smooth1), 3)

CP2 <- mean(cover2)
SCP2 <- mean(cover_smooth2)
SE2 <- round(mean(Std2), 3)
SSE2 <- round(mean(Std_smooth2), 3)

CR <- round(mean(censor_rate), 3)
V_test_opt <- round(V_test_opt, 3)


colNames <- c('CR', 'PS', 'eta0_hat', 'eta1_hat', 'eta2_hat', 'MRL_hat(eta_hat)', 'SE', 'SE1', 'SE2', 'CP', 'CP1', 'CP2', 'MRL(eta_hat)', 'CA', 'MRL', 'sd(Mhat_opt)', 'sd(Mhat_opt1)')
rowNames <- c('IPWE', 'SIPWE')
row1 <- c(CR, PS, eta0_hat, eta1_hat, eta2_hat, MRL_hat, SE, SE1, SE2, CP, CP1, CP2, MRL, CA, V_test_opt, sd(Mhat_opt), sd(Mhat_opt1))
row2 <- c(CR, PS, Seta0_hat, Seta1_hat, Seta2_hat, SMRL_hat, SSE, SSE1, SSE2, SCP, SCP1, SCP2, SMRL, SCA, V_test_opt, sd(Mhat_opt_smooth), sd(Mhat_opt_smooth1))
results <- data.frame(rbind(row1, row2), row.names = rowNames)
colnames(results) <- colNames

fileName <- paste('time-varying regime', 'IPW', 'proportional', N, PS, t, CR, sep = '-')
fileName <- paste(fileName, 'csv', sep = '.')
#fileName <- paste('/Users/liuzhishuai/Desktop/DTR-MRL/数值结果202105/', fileName, sep = '')
write.table(results, file = fileName, sep = ',', row.names = TRUE, col.names = NA)




