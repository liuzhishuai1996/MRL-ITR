rm(list=ls())
setwd()
library('rgenoud')
library('survival')
library('evd')
expit <- function(x) exp(x)/(1 + exp(x))
library(compiler)
enableJIT(3)
library(survival)
library(rgenoud)
source("RealData-Functions-v3.R")
Frd <- function(t0, DataList, Neta, PS)
{ 
  if (PS == TRUE){
    ps <- Fps.true(DataList)
  }else{
    ps <- Fps(DataList) 
  }
  reg <- Freg(DataList)
  PC <- FPC(DataList)
  prep <- Fprep(DataList, t0, ps, reg, PC)
  
  result <- matrix(NA, nrow=4, ncol=3, dimnames=list(c("I","A", "SI", "SA"), c("int", "X1", "X2")))
  
  result[1, ] <- Genetic.IPWE(t0, DataList, ps, TRUE, Neta)[1:Neta]
  
  result[2, ] <- Genetic.AIPWE(t0, DataList, prep, TRUE, Neta)[1:Neta]
  
  result[3, ] <- Genetic.IPWE(t0, DataList, ps, FALSE, Neta)[1:Neta]
  
  result[4, ] <- Genetic.AIPWE(t0, DataList, prep, FALSE, Neta)[1:Neta]
  
  return(result)
}

N <- 500
t <- 2
p <- 1; q <- - 3; c <- 0.5; d <- 1
C0 <- 0 
C1 <- 17 
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
  
  
  data1 <- data.frame(delta = as.integer(delta), Ttilde = Y, A = treat, X1 = X1, X2 = X2)
  index <- order(data1$Ttilde, decreasing=F)  
  data <- data1[index,]
  
  DataList <- list()
  DataList$Ttilde <- data$Ttilde
  DataList$delta <- data$delta
  DataList$A <- data$A
  DataList$X <- model.matrix( ~ X1+X2, data)[,-1]
  Neta <- ncol(DataList$X)+1
  
  
  
  pre0 <- matrix(NA, nrow=4, ncol=3, 
                 dimnames=list(paste0(rep(c(t), each=4), ",", c("I","A","SI","SA")),
                               c("int", "X1", "X2")))
  
  pre0[1:4, ] <- Frd(t0=t, DataList, Neta, correctPSmodel)
  
  ETA_IPW[iter, ] <- pre0[1,]
  ETA_AIPW[iter, ] <- pre0[2,]
  SETA_IPW[iter, ] <- pre0[3,]
  SETA_AIPW[iter, ] <- pre0[4,]
  
  treat_IPW_test <- as.numeric(as.matrix(cbind(rep(1, N_test), Xclassify_test)) %*% pre0[1,] > 0)
  treat_AIPW_test <- as.numeric(as.matrix(cbind(rep(1, N_test), Xclassify_test)) %*% pre0[2,] > 0)
  Streat_IPW_test <- as.numeric(as.matrix(cbind(rep(1, N_test), Xclassify_test)) %*% pre0[3,] > 0)
  Streat_AIPW_test <- as.numeric(as.matrix(cbind(rep(1, N_test), Xclassify_test)) %*% pre0[4,] > 0)
  
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
  V_IPW_test[iter] <- sum(T_IPW_test[T_IPW_test > t] - t)/sum(T_IPW_test > t)
  V_AIPW_test[iter] <- sum(T_AIPW_test[T_AIPW_test > t] - t)/sum(T_AIPW_test > t)
  SV_IPW_test[iter] <- sum(ST_IPW_test[ST_IPW_test > t] - t)/sum(ST_IPW_test > t)
  SV_AIPW_test[iter] <- sum(ST_AIPW_test[ST_AIPW_test > t] - t)/sum(ST_AIPW_test > t)
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

fileName <- paste('time-varying regime', 'SP', N, PS, t, CR, sep = '-')
fileName <- paste(fileName, 'csv', sep = '.')
#fileName <- paste('/Users/liuzhishuai/Documents/硕士/研三上学期/DTR-MRL/Review/data_codes_survival/', fileName, sep = '')
write.table(results, file = fileName, sep = ',', row.names = TRUE, col.names = NA)
