# Basic functions for Simulation and Real data analysis
# Zhishuai Liu
MRL <- function(data, TIME, N2, D1, D2, type){
  #========proportional=========
  if (type == "proportional"){
    index <- order(data$Y)
    data <- data[index, ]
    
    #==========step 0: Initialization========#
    loop <- 0
    Bt <- 0
    INVBt <- 0
    CRITERION <- 0.0000000001
    
    Exx <- matrix(0, nrow = N2 + 1, ncol = 1)
    Ezx <- matrix(0, nrow = N2 + 1, ncol = 3)
    Ezz <- matrix(0, nrow = N2 + 1, ncol = 9)
    
    tempxx <- NA; tempzx <- rep(NA, 3); tempzz <- rep(NA, 9)
    temp <- NA; temp1 <- NA; tempbeta <- rep(NA, 3); tempalpha <- NA; A <- matrix(NA, nrow = 3, ncol = 3)
    tempest <- rep(0, 4)
    tempest[1] <- log(D1*TIME + D2); tempest[2] <- -0.5; tempest[3] <- 0.5; tempest[4] <- -0.5
    oldtempest <- rep(0, 4)
    oldalpha <- rep(0, N2 + 1) 
    alpha <- log(D1*c(0, data$Y) + D2)
    diff <- 100
    
    #==========step 1: Find the position t of the last invertible Exx========#
    
    # indi <- 1
    # t <- 0
    # while(indi > 0.5 && t <= (N2 - 1)){
    #   tempxx <- 0
    #   
    #   for (i in (t + 1):N2){
    #     temp <- alpha[t + 1] + tempest[2] * data$Z1[i] + tempest[3] * data$AZ1[i] + tempest[4] * data$AZ2[i]
    #     temp <- data$delta[i]*exp(temp)/data$Sc[i]/N2
    #     
    #     tempxx <- tempxx + temp
    #   }
    #   
    #   if (abs(tempxx) < CRITERION){
    #     indi <- 0
    #     Singular <- t - 1
    #   }
    #   
    #   t <- t + 1
    # }
    
    indi <- 2
    t <- N2
    while(indi > 0.5 && t >= 1){
      if (data$delta[t] > 0.5){
        indi <- indi - 1
      }
      t <- t - 1
    }
    indi <- t - 5 
    
    #==========step 2: Iteration for beta and alpha========#
    while(diff > 0.000001 & loop < 100){
      loop <- loop + 1
      
      cat('The', loop,'th iteration', '\n')
      
      # 2.1 Re-initialization, that is, make "new" "old".
      for (t in 0:indi){
        oldalpha[t + 1] <- alpha[t + 1]
      }
      
      oldtempest <- tempest
      
      tempxx <- 0; tempzx <- rep(0, 3); tempzz <- rep(0, 9)
      temp <- 0
      
      # 2.2: Compute Exx, Ezx and Ezz. And then new beta and A.
      
      tempbeta <- rep(0, 3)
      A <- matrix(0, nrow = 3, ncol = 3)
      
      #Integration over H(t) for t = 1 to indi.
      for (t in 1:indi){
        # Compute Exx, Ezx and Ezz
        
        tempxx <- 0; tempzx <- rep(0, 3); tempzz <- rep(0, 9)
        temp <- 0
        
        # for (i in (t + 1):N2){
        #   temp <- data$delta[i] * exp(oldalpha[t + 1] + oldtempest[2] * data$Z1[i] + oldtempest[3] * data$AZ1[i] + oldtempest[4] * data$AZ2[i]) * data$Sc[i] / N2
        #   
        #   tempxx <- tempxx + temp
        #   
        #   tempzx <- tempzx + temp*c(data$Z1[i], data$AZ1[i], data$AZ2[i])
        #   
        #   tempzz <- tempzz + temp*as.vector(c(data$Z1[i], data$AZ1[i], data$AZ2[i]) %*% t(c(data$Z1[i], data$AZ1[i], data$AZ2[i])))
        # }
        temp <- data$delta[(t + 1):N2]*exp(oldalpha[t + 1] + oldtempest[2] * data$Z1[(t + 1):N2] + oldtempest[3] * data$AZ1[(t + 1):N2] + oldtempest[4] * data$AZ2[(t + 1):N2])*data$Sc[(t + 1):N2]/N2
        tempxx <- sum(temp)
        tempzx[1] <- sum(temp*data$Z1[(t + 1):N2])
        tempzx[2] <- sum(temp*data$AZ1[(t + 1):N2])
        tempzx[3] <- sum(temp*data$AZ2[(t + 1):N2])
        tempzz <- matrix(t(temp*cbind(data$Z1[(t + 1):N2], data$AZ1[(t + 1):N2], data$AZ2[(t + 1):N2])) %*% cbind(data$Z1[(t + 1):N2], data$AZ1[(t + 1):N2], data$AZ2[(t + 1):N2]), nrow = 3, ncol = 3)
        
        
        Ezz[t + 1, ] <- tempzz
        Ezx[t + 1, ] <- tempzx
        Exx[t + 1] <- tempxx 
        
        Bt <- tempxx
        INVBt <- 1/Bt
        
        A <- A + matrix(Ezz[t + 1, ], nrow = 3, ncol = 3) - INVBt*(Ezx[t + 1, ] %*% t(Ezx[t + 1, ]))
        
        # for (i in (t + 1):N2){
        #   temp <- data$delta[i]*(data$Y[i] - data$Y[t] - exp(oldalpha[t + 1] + oldtempest[2] * data$Z1[i] + oldtempest[3] * data$AZ1[i] + oldtempest[4] * data$AZ2[i]))*data$Sc[i]
        #   temp <- temp*(c(data$Z1[i], data$AZ1[i], data$AZ2[i]) - Ezx[t + 1, ]*INVBt*data$X[i])
        #   
        #   tempbeta <- tempbeta + temp
        # }
        temp <- data$delta[(t + 1):N2]*(data$Y[(t + 1):N2] - data$Y[t] - exp(oldalpha[t + 1] + oldtempest[2] * data$Z1[(t + 1):N2] + oldtempest[3] * data$AZ1[(t + 1):N2] + oldtempest[4] * data$AZ2[(t + 1):N2]))*data$Sc[(t + 1):N2]
        tempbeta[1] <- sum(temp*(data$Z1[(t + 1):N2] - Ezx[t + 1, 1]*INVBt*data$X[(t + 1):N2])) + tempbeta[1]
        tempbeta[2] <- sum(temp*(data$AZ1[(t + 1):N2] - Ezx[t + 1, 2]*INVBt*data$X[(t + 1):N2])) + tempbeta[2]
        tempbeta[3] <- sum(temp*(data$AZ2[(t + 1):N2] - Ezx[t + 1, 3]*INVBt*data$X[(t + 1):N2])) + tempbeta[3]
      }
      
      tempest[2:4] <- oldtempest[2:4] + tempbeta %*% solve(A) / N2
      
      # 2.3: Compute alpha[t][] for all t (both censored and uncensored, from 0 to indi).
      for (t in 1:indi){
        Bt <- Exx[t + 1]
        
        INVBt <- 1/Bt
        
        tempalpha <- 0
        
        # for (i in (t + 1):N2){
        #   temp <- oldalpha[t + 1] + oldtempest[2] * data$Z1[i] + oldtempest[3] * data$AZ1[i] + oldtempest[4] * data$AZ2[i]
        #   temp <- exp(temp) + exp(temp)*t(c(data$Z1[i], data$AZ1[i], data$AZ2[i])) %*% (tempest[2:4] - oldtempest[2:4])
        #   temp <- data$delta[i]*(data$Y[i] - data$Y[t] - temp)*data$Sc[i]/N2
        #   tempalpha <- tempalpha + temp
        # }
        temp <- exp(oldalpha[t + 1] + oldtempest[2] * data$Z1[(t + 1):N2] + oldtempest[3] * data$AZ1[(t + 1):N2] + oldtempest[4] * data$AZ2[(t + 1):N2])
        aaa <- temp + temp*data$Z1[(t + 1):N2]*(tempest[2] - oldtempest[2]) + temp*data$AZ1[(t + 1):N2]*(tempest[3] - oldtempest[3]) + temp*data$AZ2[(t + 1):N2]*(tempest[4] - oldtempest[4])
        bbb <- data$delta[(t + 1):N2]*(data$Y[(t + 1):N2] - data$Y[t] - aaa)*data$Sc[(t + 1):N2]/N2
        tempalpha <- sum(bbb)
        
        
        alpha[t + 1] <- oldalpha[t + 1] + INVBt*tempalpha
      }
      
      # 2.4: Compute tempest
      
      tempxx <- 0
      temp <- temp1 <- 0
      tempalpha <- 0
      
      # for (i in 1:N2){
      #   if (data$Y[i] > TIME){
      #     temp <- oldtempest[1] + oldtempest[2] * data$Z1[i] + oldtempest[3] * data$AZ1[i] + oldtempest[4] * data$AZ2[i]
      #     
      #     temp1 <- exp(temp) + exp(temp)*t(c(data$Z1[i], data$AZ1[i], data$AZ2[i])) %*% (tempest[2:4] - oldtempest[2:4])
      #     temp1 <- data$delta[i]*(data$Y[i] - TIME - temp1)*data$Sc[i]/N2
      #     tempalpha <- tempalpha + temp1
      #     
      #     temp1 <- data$delta[i]*exp(temp)*data$Sc[i]/N2
      #     tempxx <- tempxx + temp1
      #   }
      # }
      temp <- exp(oldtempest[1] + oldtempest[2] * data$Z1[1:N2] + oldtempest[3] * data$AZ1[1:N2] + oldtempest[4] * data$AZ2[1:N2])
      aaa <- temp + temp*data$Z1[1:N2]*(tempest[2] - oldtempest[2]) + temp*data$AZ1[1:N2]*(tempest[3] - oldtempest[3]) + temp*data$AZ2[1:N2]*(tempest[4] - oldtempest[4])
      bbb <- as.numeric(data$Y[1:N2] > TIME)*data$delta[1:N2]*(data$Y[1:N2] - TIME - aaa)*data$Sc[1:N2]/N2
      tempalpha <- sum(bbb)
      
      ccc <- as.numeric(data$Y[1:N2] > TIME)*data$delta[1:N2]*temp*data$Sc[1:N2]/N2
      tempxx <- sum(ccc)
      
      
      if (abs(tempxx) > CRITERION){
        Bt <- tempxx
        INVBt <- 1/Bt
        
        tempest[1] = oldtempest[1] + INVBt*tempalpha
      }else{
        print("This t is not estimable since the matrix is singular")
      }
      
      diff <- 0
      
      for (t in 0:3){
        diff <- diff + abs(tempest[t + 1] - oldtempest[t + 1])
      }
      
    }
  }
  #========additive=========
  if (type == "additive"){
    index <- order(data$Y)
    data <- data[index, ]
    
    #==========step 0: Initialization========#
    loop <- 0
    Bt <- 0
    INVBt <- 0
    CRITERION <- 0.0000000001
    
    Exx <- matrix(0, nrow = N2 + 1, ncol = 1)
    Ezx <- matrix(0, nrow = N2 + 1, ncol = 3)
    Ezz <- matrix(0, nrow = N2 + 1, ncol = 9)
    
    tempxx <- NA; tempzx <- rep(NA, 3); tempzz <- rep(NA, 9)
    temp <- NA; temp1 <- NA; tempbeta <- rep(NA, 3); tempalpha <- NA; A <- matrix(NA, nrow = 3, ncol = 3)
    tempest <- rep(0, 4)
    tempest[1] <- D1*TIME + D2; tempest[2] <- - 0.5; tempest[3] <- 0.5; tempest[4] <- - 0.5
    oldtempest <- rep(0, 4)
    oldalpha <- rep(0, N2 + 1) 
    # alpha <- rep(0, N2 + 1) 
    alpha <- D1*c(0, data$Y) + D2
    diff <- 100
    
    #==========step 1: Find the position t of the last invertible Exx========#
    
    # indi <- 1
    # t <- 0
    # while(indi > 0.5 && t <= (N2 - 1)){
    #   tempxx <- 0
    #   
    #   for (i in (t + 1):N2){
    #     temp <- alpha[t + 1] + tempest[2] * data$Z1[i] + tempest[3] * data$AZ1[i] + tempest[4] * data$AZ2[i]
    #     temp <- data$delta[i]/data$Sc[i]/N2
    #     
    #     tempxx <- tempxx + temp
    #   }
    #   
    #   if (abs(tempxx) < CRITERION){
    #     indi <- 0
    #     Singular <- t - 1
    #   }
    #   
    #   t <- t + 1
    # }
    
    indi <- 2
    t <- N2
    while(indi > 0.5 && t >= 1){
      if (data$delta[t] > 0.5){
        indi <- indi - 1
      }
      t <- t - 1
    }
    indi <- t - 5 
    
    #==========step 2: Iteration for beta and alpha========#
    while(diff > 0.0001 & loop < 100){
      loop <- loop + 1
      
      cat('The', loop,'th iteration', '\n')
      
      # 2.1 Re-initialization, that is, make "new" "old".
      for (t in 0:indi){
        oldalpha[t + 1] <- alpha[t + 1]
      }
      
      oldtempest <- tempest
      
      tempxx <- 0; tempzx <- rep(0, 3); tempzz <- rep(0, 9)
      temp <- 0
      
      # 2.2: Compute Exx, Ezx and Ezz. And then new beta and A.
      
      tempbeta <- rep(0, 3)
      A <- matrix(0, nrow = 3, ncol = 3)
      
      #Integration over H(t) for t = 1 to indi.
      for (t in 1:indi){ 
        # Compute Exx, Ezx and Ezz
        
        tempxx <- 0; tempzx <- rep(0, 3); tempzz <- rep(0, 9)
        temp <- 0
        
        # for (i in (t + 1):N2){
        #   temp <- oldalpha[t + 1] + oldtempest[2] * data$Z1[i] + oldtempest[3] * data$AZ1[i] + oldtempest[4] * data$AZ2[i]
        #   temp <- data$delta[i] * temp * data$Sc[i] / N2
        #   
        #   tempxx <- tempxx + temp
        #   
        #   tempzx <- tempzx + temp*c(data$Z1[i], data$AZ1[i], data$AZ2[i])
        #   
        #   tempzz <- tempzz + temp*as.vector(c(data$Z1[i], data$AZ1[i], data$AZ2[i]) %*% t(c(data$Z1[i], data$AZ1[i], data$AZ2[i])))
        # }
        temp <- data$delta[(t + 1):N2]*(oldalpha[t + 1] + oldtempest[2] * data$Z1[(t + 1):N2] + oldtempest[3] * data$AZ1[(t + 1):N2] + oldtempest[4] * data$AZ2[(t + 1):N2])*data$Sc[(t + 1):N2]/N2
        tempxx <- sum(temp)
        tempzx[1] <- sum(temp*data$Z1[(t + 1):N2])
        tempzx[2] <- sum(temp*data$AZ1[(t + 1):N2])
        tempzx[3] <- sum(temp*data$AZ2[(t + 1):N2])
        tempzz <- matrix(t(temp*cbind(data$Z1[(t + 1):N2], data$AZ1[(t + 1):N2], data$AZ2[(t + 1):N2])) %*% cbind(data$Z1[(t + 1):N2], data$AZ1[(t + 1):N2], data$AZ2[(t + 1):N2]), nrow = 3, ncol = 3)
        
        Ezz[t + 1, ] <- tempzz
        Ezx[t + 1, ] <- tempzx
        Exx[t + 1] <- tempxx 
        
        Bt <- tempxx
        INVBt <- 1/Bt
        
        A <- A + matrix(Ezz[t + 1, ], nrow = 3, ncol = 3) - INVBt*(Ezx[t + 1, ] %*% t(Ezx[t + 1, ]))
        
        # for (i in (t + 1):N2){
        #   temp <- oldalpha[t + 1] + oldtempest[2] * data$Z1[i] + oldtempest[3] * data$AZ1[i] + oldtempest[4] * data$AZ2[i]
        #   temp <- data$delta[i]*(data$Y[i] - data$Y[t] - temp)*data$Sc[i]
        #   temp <- temp*(c(data$Z1[i], data$AZ1[i], data$AZ2[i]) - Ezx[t + 1, ]*INVBt*data$X[i])
        #   
        #   tempbeta <- tempbeta + temp
        # }
        temp <- data$delta[(t + 1):N2]*(data$Y[(t + 1):N2] - data$Y[t] - (oldalpha[t + 1] + oldtempest[2] * data$Z1[(t + 1):N2] + oldtempest[3] * data$AZ1[(t + 1):N2] + oldtempest[4] * data$AZ2[(t + 1):N2]))*data$Sc[(t + 1):N2]
        tempbeta[1] <- sum(temp*(data$Z1[(t + 1):N2] - Ezx[t + 1, 1]*INVBt*data$X[(t + 1):N2])) + tempbeta[1]
        tempbeta[2] <- sum(temp*(data$AZ1[(t + 1):N2] - Ezx[t + 1, 2]*INVBt*data$X[(t + 1):N2])) + tempbeta[2]
        tempbeta[3] <- sum(temp*(data$AZ2[(t + 1):N2] - Ezx[t + 1, 3]*INVBt*data$X[(t + 1):N2])) + tempbeta[3]
      }
      
      tempest[2:4] <- oldtempest[2:4] + tempbeta %*% solve(A) / N2
      
      # 2.3: Compute alpha[t][] for all t (both censored and uncensored, from 0 to indi).
      for (t in 1:indi){
        Bt <- Exx[t + 1]
        
        INVBt <- 1/Bt
        
        tempalpha <- 0
        
        # for (i in (t + 1):N2){
        #   temp <- oldalpha[t + 1] + oldtempest[2] * data$Z1[i] + oldtempest[3] * data$AZ1[i] + oldtempest[4] * data$AZ2[i]
        #   temp <- temp + t(c(data$Z1[i], data$AZ1[i], data$AZ2[i])) %*% (tempest[2:4] - oldtempest[2:4])
        #   temp <- data$delta[i]*(data$Y[i] - data$Y[t] - temp)*data$Sc[i]/N2
        #   tempalpha <- tempalpha + temp
        # }
        temp <- oldalpha[t + 1] + oldtempest[2] * data$Z1[(t + 1):N2] + oldtempest[3] * data$AZ1[(t + 1):N2] + oldtempest[4] * data$AZ2[(t + 1):N2]
        aaa <- temp + data$Z1[(t + 1):N2]*(tempest[2] - oldtempest[2]) + data$AZ1[(t + 1):N2]*(tempest[3] - oldtempest[3]) + data$AZ2[(t + 1):N2]*(tempest[4] - oldtempest[4])
        bbb <- data$delta[(t + 1):N2]*(data$Y[(t + 1):N2] - data$Y[t] - aaa)*data$Sc[(t + 1):N2]/N2
        tempalpha <- sum(bbb)
        
        alpha[t + 1] <- oldalpha[t + 1] + INVBt*tempalpha
      }
      
      # 2.4: Compute tempest
      
      tempxx <- 0
      temp <- temp1 <- 0
      tempalpha <- 0
      
      # for (i in 1:N2){
      #   if (data$Y[i] > TIME){
      #     temp <- oldtempest[1] + oldtempest[2] * data$Z1[i] + oldtempest[3] * data$AZ1[i] + oldtempest[4] * data$AZ2[i]
      #     
      #     temp1 <- temp + t(c(data$Z1[i], data$AZ1[i], data$AZ2[i])) %*% (tempest[2:4] - oldtempest[2:4])
      #     temp1 <- data$delta[i]*(data$Y[i] - TIME - temp1)*data$Sc[i]/N2
      #     tempalpha <- tempalpha + temp1
      #     
      #     temp1 <- data$delta[i]*data$Sc[i]/N2
      #     tempxx <- tempxx + temp1
      #   }
      # }
      temp <- (oldtempest[1] + oldtempest[2] * data$Z1[1:N2] + oldtempest[3] * data$AZ1[1:N2] + oldtempest[4] * data$AZ2[1:N2])
      aaa <- temp + data$Z1[1:N2]*(tempest[2] - oldtempest[2]) + data$AZ1[1:N2]*(tempest[3] - oldtempest[3]) + data$AZ2[1:N2]*(tempest[4] - oldtempest[4])
      bbb <- as.numeric(data$Y[1:N2] > TIME)*data$delta[1:N2]*(data$Y[1:N2] - TIME - aaa)*data$Sc[1:N2]/N2
      tempalpha <- sum(bbb)
      
      ccc <- as.numeric(data$Y[1:N2] > TIME)*data$delta[1:N2]*temp*data$Sc[1:N2]/N2
      tempxx <- sum(ccc)
      
      if (abs(tempxx) > CRITERION){
        Bt <- tempxx
        INVBt <- 1/Bt
        
        tempest[1] = oldtempest[1] + INVBt*tempalpha
      }else{
        print("This t is not estimable since the matrix is singular")
      }
      
      diff <- 0
      
      for (t in 0:3){
        diff <- diff + abs(tempest[t + 1] - oldtempest[t + 1])
      }
      
    }
  }
  #=======return the estimator=======
  result <- list()
  result$est <- tempest
  result$alpha <- alpha
  result$indi <- indi
  return(result)
}

MRL_ESTSTD <- function(data, alpha, TIME_seq, N2, phihat, type, tempest){
  index <- order(data$Y)
  data <- data[index, ]
  
  if (type == "proportional"){
    
    # Step 0: initialization
    
    Exx <- matrix(0, nrow = N2 + 1, ncol = 1)
    Ezx <- matrix(0, nrow = N2 + 1, ncol = 3)
    Ezz <- matrix(0, nrow = N2 + 1, ncol = 9)
    
    
    mit <- matrix(0, nrow = N2 + 1, ncol = N2 + 1)
    mzeh <- matrix(0, nrow = N2 + 1, ncol = 3); eta <- matrix(0, nrow = N2 + 1, ncol = 3); xi <- matrix(0, nrow = N2 + 1, ncol = 3)
    tempExx <- 0; tempExxInv <- 0
    
    alpha_mit <- rep(0, N2 + 1); alpha_q <- rep(0, N2 + 1); alpha_Exx <- 0; alpha_Ezx <- rep(0, 3)
    phi_alpha <- matrix(0, nrow = N2 + 1, ncol = length(TIME_seq))
    # Step 1: Compute Ezz, Ezx, Exx
    
    for (t in 1:indi){
      tempxx <- 0; tempzx <- rep(0, 3); tempzz <- rep(0, 9)
      temp <- 0
      # 
      # for (i in (t + 1):N2){
      #   temp <- exp(alpha[t + 1] + tempest[2] * data$Z1[i] + tempest[3] * data$AZ1[i] + tempest[4] * data$AZ2[i])
      #   
      #   tempxx <- tempxx + temp*data$delta[i]*data$Sc[i]/N2
      #   
      #   tempzx <- tempzx + temp*data$delta[i]*c(data$Z1[i], data$AZ1[i], data$AZ2[i])
      #   
      #   tempzz <- tempzz + temp*data$delta[i]*as.vector(c(data$Z1[i], data$AZ1[i], data$AZ2[i]) %*% t(c(data$Z1[i], data$AZ1[i], data$AZ2[i])))
      # }
      temp <- exp(alpha[t + 1] + tempest[2] * data$Z1[(t + 1):N2] + tempest[3] * data$AZ1[(t + 1):N2] + tempest[4] * data$AZ2[(t + 1):N2])
      tempxx <- sum(temp*data$delta[(t + 1):N2]*data$Sc[(t + 1):N2]/N2)
      tempzx[1] <- sum(temp*data$delta[(t + 1):N2]*data$Z1[(t + 1):N2]/N2)
      tempzx[2] <- sum(temp*data$delta[(t + 1):N2]*data$AZ1[(t + 1):N2]/N2)
      tempzx[3] <- sum(temp*data$delta[(t + 1):N2]*data$AZ2[(t + 1):N2]/N2)
      tempzz <- matrix(t(temp*data$delta[(t + 1):N2]*cbind(data$Z1[(t + 1):N2], data$AZ1[(t + 1):N2], data$AZ2[(t + 1):N2])) %*% cbind(data$Z1[(t + 1):N2], data$AZ1[(t + 1):N2], data$AZ2[(t + 1):N2]), nrow = 3, ncol = 3)/N2
      
      Ezz[t + 1, ] <- tempzz
      Ezx[t + 1, ] <- tempzx
      Exx[t + 1] <- tempxx
    }
    
    # Step 2: Compute mit
    for (i in 1:N2){ # i > t
      
      # temp <- 0
      # 
      # for (t in 1:(i - 1)){ # no need to worry about t=0 since H(t) doesn't jump at t=0.
      #   
      #   if (t <= indi){
      #     
      #     temp <- exp(alpha[t + 1] + tempest[2] * data$Z1[i] + tempest[3] * data$AZ1[i] + tempest[4] * data$AZ2[i])
      #     mit[i, t + 1] <- data$delta[i]*data$Sc[i]*(data$Y[i] - data$Y[t + 1] - temp)
      #     
      #   }
      # }
      temp <- exp(alpha[2:i] + tempest[2] * data$Z1[i] + tempest[3] * data$AZ1[i] + tempest[4] * data$AZ2[i])
      mit[i, 2:i] <- (t <= indi)*data$delta[i]*data$Sc[i]*(data$Y[i] - data$Y[2:i] - temp)
    }
    
    # Step 3: Compute eta, Sigma
    
    # First we need to compute mzeh -- the first part of eta. This part does not depend on v and can be
    # computed and stored before computing eta to save time. It depends only on i.
    for (i in 2:N2){
      # temp <- rep(0, 3)
      # for (t in 1:(i - 1)){# if t >= i then mit is 0 so no need to proceed.
      #   if (t <= indi & data$delta[t] > 0.5){
      #     tempExx <- Exx[t + 1]
      #     tempExxInv <- 1/tempExx
      #     tempzx <- Ezx[t + 1, ]*tempExxInv
      #     temp <- temp + mit[i, t + 1]*(c(data$Z1[i], data$AZ1[i], data$AZ2[i]) - tempzx)
      #   }
      # }
      
      temp <- ( as.numeric((data$delta[1:(min(indi + 1, i) - 1)])) * mit[i, 2:min(indi + 1, i)] ) %*% ( matrix(rep(c(data$Z1[i], data$AZ1[i], data$AZ2[i]), min(indi + 1, i) - 1), nrow = min(indi + 1, i) - 1, ncol = 3, byrow = T) - 1/Exx[2:min(indi + 1, i)]*Ezx[2:min(indi + 1, i), ])
      
      # mzeh1[i, ] <- temp
      mzeh[i, ] <- temp
      
    }
    
    # Second for eta[V_j], where V_j is indexed by j
    for (j in 1:N2){
      
      temp <- rep(0, 3)
      # for (i in 1:N2){
      #   temp <- temp + mzeh[i, ]*phihat[j, i]*data$Sc[i] 	# for phihat: the j-th obs, time is i
      # }
      temp[1] <- sum(mzeh[1:N2, 1]*phihat[j, ]*data$Sc)
      temp[2] <- sum(mzeh[1:N2, 2]*phihat[j, ]*data$Sc)
      temp[3] <- sum(mzeh[1:N2, 3]*phihat[j, ]*data$Sc)
      
      
      eta[j, ] <- temp/N2
      
      # now compute xi[j]
      xi[j, ] = mzeh[j, ] - eta[j, ];
    }
    
    # Step 4:	compute A and the std of beta
    
    # A <- matrix(0, nrow = 3, ncol = 3)
    # for (t in 1:indi){
    #   if (data$delta[t] > 0.5){ # if this t is uncensored.
    #     tempExx <- Exx[t + 1]
    # 
    #     tempExxInv <- 1/tempExx
    # 
    #     A <- A + matrix(Ezz[t + 1, ], nrow = 3, ncol = 3) - tempExxInv*(Ezx[t + 1, ] %*% t(Ezx[t + 1, ]))
    #   }
    # }
    A <- matrix(apply(data$delta[1:indi]*Ezz[2:(indi + 1), ], 2, sum), nrow = 3, ncol = 3) - t(data$delta[1:indi]/Exx[2:(indi + 1)]*Ezx[2:(indi + 1), ]) %*% (Ezx[2:(indi + 1), ])
    
    
    phi_beta <- t(solve(A) %*% t(xi))
    
    # Step 5: Compute the std of alpha[t] for t = TIME
    for (o in 1:length(TIME_seq)){
      count  <- 1 #count is used to index the location of TIME in data
      
      while (data$Y[count] <= TIME_seq[o] & count < N2){
        count <- count + 1
      } 
      # So now count is the position with which Exx, Ezx and mit should start.
      
      # First compute alpha.mit[i]
      # for (i in count:N2){
      #   temp <- exp(tempest[1] + tempest[2] * data$Z1[i] + tempest[3] * data$AZ1[i] + tempest[4] * data$AZ2[i])
      #   alpha_mit[i] <- data$delta[i]*data$Sc[i]*(data$Y[i] - TIME - temp)
      # }
      temp <- exp(tempest[1] + tempest[2] * data$Z1[count:N2] + tempest[3] * data$AZ1[count:N2] + tempest[4] * data$AZ2[count:N2])
      alpha_mit[count:N2] <- data$delta[count:N2]*data$Sc[count:N2]*(data$Y[count:N2] - TIME_seq[o] - temp)
      
      # Second compute alpha.Exx and alpha.Exz
      tempxx <- 0; tempzx <- rep(0, 3)
      temp <- 0
      
      # for (i in count:N2){ # Y_i > t
      #   temp <- exp(tempest[1] + tempest[2] * data$Z1[i] + tempest[3] * data$AZ1[i] + tempest[4] * data$AZ2[i])
      #   
      #   tempxx <- tempxx + temp*data$delta[i]*data$Sc[i]
      #   
      #   tempzx <- tempzx + temp*data$delta[i]*data$Sc[i]*c(data$Z1[i], data$AZ1[i], data$AZ2[i])
      # }
      temp <- exp(tempest[1] + tempest[2] * data$Z1[count:N2] + tempest[3] * data$AZ1[count:N2] + tempest[4] * data$AZ2[count:N2])
      tempxx <- sum(temp*data$delta[count:N2]*data$Sc[count:N2])
      tempzx[1] <- sum(temp*data$delta[count:N2]*data$Sc[count:N2]*data$Z1[count:N2])
      tempzx[2] <- sum(temp*data$delta[count:N2]*data$Sc[count:N2]*data$AZ1[count:N2])
      tempzx[3] <- sum(temp*data$delta[count:N2]*data$Sc[count:N2]*data$AZ2[count:N2])
      
      alpha_Exx <- tempxx/N2
      
      alpha_Ezx <- tempzx/N2
      
      # Third compute alpha.q[j]
      tempExx <- alpha_Exx
      tempExxInv <- 1/tempExx
      for (j in 1:N2){
        tempxx <- 0
        # for (i in 1:N2){
        #   tempxx <- tempxx + alpha_mit[i]*phihat[j, i]*data$Sc[i]
        # }
        tempxx <- sum(alpha_mit[1:N2]*phihat[j, ]*data$Sc)
        alpha_q[j] <- tempxx/N2
        
        # Now everything is ready.
        
        tempxx <- alpha_mit[j] - alpha_q[j] - alpha_Ezx %*% t(xi[j, ] %*% solve(A))
        
        phi_alpha[j, o] <- tempExxInv*tempxx
      }
    }
    
  }
  
  if (type == "additive"){
    # Step 0: initialization
    
    Exx <- matrix(0, nrow = N2 + 1, ncol = 1)
    Ezx <- matrix(0, nrow = N2 + 1, ncol = 3)
    Ezz <- matrix(0, nrow = N2 + 1, ncol = 9)
    
    
    mit <- matrix(0, nrow = N2 + 1, ncol = N2 + 1)
    mzeh <- matrix(0, nrow = N2 + 1, ncol = 3); eta <- matrix(0, nrow = N2 + 1, ncol = 3); xi <- matrix(0, nrow = N2 + 1, ncol = 3)
    tempExx <- 0; tempExxInv <- 0
    
    alpha_mit <- rep(0, N2 + 1); alpha_q <- rep(0, N2 + 1); alpha_Exx <- 0; alpha_Ezx <- rep(0, 3)
    phi_alpha <- matrix(0, nrow = N2 + 1, ncol = length(TIME_seq))
    # Step 1: Compute Ezz, Ezx, Exx
    
    for (t in 1:indi){
      tempxx <- 0; tempzx <- rep(0, 3); tempzz <- rep(0, 9)
      temp <- 0
      # 
      # for (i in (t + 1):N2){
      #   temp <- exp(alpha[t + 1] + tempest[2] * data$Z1[i] + tempest[3] * data$AZ1[i] + tempest[4] * data$AZ2[i])
      #   
      #   tempxx <- tempxx + temp*data$delta[i]*data$Sc[i]/N2
      #   
      #   tempzx <- tempzx + temp*data$delta[i]*c(data$Z1[i], data$AZ1[i], data$AZ2[i])
      #   
      #   tempzz <- tempzz + temp*data$delta[i]*as.vector(c(data$Z1[i], data$AZ1[i], data$AZ2[i]) %*% t(c(data$Z1[i], data$AZ1[i], data$AZ2[i])))
      # }
      
      # temp <- exp(alpha[t + 1] + tempest[2] * data$Z1[(t + 1):N2] + tempest[3] * data$AZ1[(t + 1):N2] + tempest[4] * data$AZ2[(t + 1):N2])
      temp <- alpha[t + 1] + tempest[2] * data$Z1[(t + 1):N2] + tempest[3] * data$AZ1[(t + 1):N2] + tempest[4] * data$AZ2[(t + 1):N2]
      tempxx <- sum(temp*data$delta[(t + 1):N2]*data$Sc[(t + 1):N2]/N2)
      tempzx[1] <- sum(temp*data$delta[(t + 1):N2]*data$Z1[(t + 1):N2])/N2
      tempzx[2] <- sum(temp*data$delta[(t + 1):N2]*data$AZ1[(t + 1):N2])/N2
      tempzx[3] <- sum(temp*data$delta[(t + 1):N2]*data$AZ2[(t + 1):N2])/N2
      tempzz <- matrix(t(temp*data$delta[(t + 1):N2]*cbind(data$Z1[(t + 1):N2], data$AZ1[(t + 1):N2], data$AZ2[(t + 1):N2])) %*% cbind(data$Z1[(t + 1):N2], data$AZ1[(t + 1):N2], data$AZ2[(t + 1):N2]), nrow = 3, ncol = 3)/N2
      
      Ezz[t + 1, ] <- tempzz
      Ezx[t + 1, ] <- tempzx
      Exx[t + 1] <- tempxx
    }
    
    # Step 2: Compute mit
    for (i in 1:N2){ # i > t
      
      # temp <- 0
      # 
      # for (t in 1:(i - 1)){ # no need to worry about t=0 since H(t) doesn't jump at t=0.
      #   
      #   if (t <= indi){
      #     
      #     temp <- exp(alpha[t + 1] + tempest[2] * data$Z1[i] + tempest[3] * data$AZ1[i] + tempest[4] * data$AZ2[i])
      #     mit[i, t + 1] <- data$delta[i]*data$Sc[i]*(data$Y[i] - data$Y[t + 1] - temp)
      #     
      #   }
      # }
      
      # temp <- exp(alpha[2:i] + tempest[2] * data$Z1[i] + tempest[3] * data$AZ1[i] + tempest[4] * data$AZ2[i])
      temp <- alpha[2:i] + tempest[2] * data$Z1[i] + tempest[3] * data$AZ1[i] + tempest[4] * data$AZ2[i] #g is the identical function 
      mit[i, 2:i] <- (t <= indi)*data$delta[i]*data$Sc[i]*(data$Y[i] - data$Y[2:i] - temp)
    }
    
    # Step 3: Compute eta, Sigma
    
    # First we need to compute mzeh -- the first part of eta. This part does not depend on v and can be
    # computed and stored before computing eta to save time. It depends only on i.
    for (i in 2:N2){
      # temp <- rep(0, 3)
      # for (t in 1:(i - 1)){# if t >= i then mit is 0 so no need to proceed.
      #   if (t <= indi & data$delta[t] > 0.5){
      #     tempExx <- Exx[t + 1]
      #     tempExxInv <- 1/tempExx
      #     tempzx <- Ezx[t + 1, ]*tempExxInv
      #     temp <- temp + mit[i, t + 1]*(c(data$Z1[i], data$AZ1[i], data$AZ2[i]) - tempzx)
      #   }
      # }
      
      temp <- ( as.numeric((data$delta[1:(min(indi + 1, i) - 1)])) * mit[i, 2:min(indi + 1, i)] ) %*% ( matrix(rep(c(data$Z1[i], data$AZ1[i], data$AZ2[i]), min(indi + 1, i) - 1), nrow = min(indi + 1, i) - 1, ncol = 3, byrow = T) - 1/Exx[2:min(indi + 1, i)]*Ezx[2:min(indi + 1, i), ])
      
      # mzeh1[i, ] <- temp
      mzeh[i, ] <- temp
      
    }
    
    # Second for eta[V_j], where V_j is indexed by j
    for (j in 1:N2){
      
      temp <- rep(0, 3)
      # for (i in 1:N2){
      #   temp <- temp + mzeh[i, ]*phihat[j, i]*data$Sc[i] 	# for phihat: the j-th obs, time is i
      # }
      # temp[1] <- sum(mzeh[1:N2, 1]*phihat[j, ]*data$Sc)
      # temp[2] <- sum(mzeh[1:N2, 2]*phihat[j, ]*data$Sc)
      # temp[3] <- sum(mzeh[1:N2, 3]*phihat[j, ]*data$Sc)
      temp[1] <- sum(mzeh[1:N2, 1]*phihat[j, ]*data$Sc)
      temp[2] <- sum(mzeh[1:N2, 2]*phihat[j, ]*data$Sc)
      temp[3] <- sum(mzeh[1:N2, 3]*phihat[j, ]*data$Sc)
      
      
      eta[j, ] <- temp/N2
      
      # now compute xi[j]
      xi[j, ] = mzeh[j, ] - eta[j, ];
    }
    
    # Step 4:	compute A and the std of beta
    
    # A <- matrix(0, nrow = 3, ncol = 3)
    # for (t in 1:indi){
    #   if (data$delta[t] > 0.5){ # if this t is uncensored.
    #     tempExx <- Exx[t + 1]
    # 
    #     tempExxInv <- 1/tempExx
    # 
    #     A <- A + matrix(Ezz[t + 1, ], nrow = 3, ncol = 3) - tempExxInv*(Ezx[t + 1, ] %*% t(Ezx[t + 1, ]))
    #   }
    # }
    A <- matrix(apply(data$delta[1:indi]*Ezz[2:(indi + 1), ], 2, sum), nrow = 3, ncol = 3) - t(data$delta[1:indi]/Exx[2:(indi + 1)]*Ezx[2:(indi + 1), ]) %*% (Ezx[2:(indi + 1), ])
    
    phi_beta <- t(solve(A) %*% t(xi))
    
    for (o in length(TIME_seq)){
      # Step 5: Compute the std of alpha[t] for t = TIME
      count  <- 1 #count is used to index the location of TIME in data
      
      while (data$Y[count] <= TIME_seq[o] & count < N2){
        count <- count + 1
      } 
      # So now count is the position with which Exx, Ezx and mit should start.
      
      # First compute alpha.mit[i]
      # for (i in count:N2){
      #   temp <- exp(tempest[1] + tempest[2] * data$Z1[i] + tempest[3] * data$AZ1[i] + tempest[4] * data$AZ2[i])
      #   alpha_mit[i] <- data$delta[i]*data$Sc[i]*(data$Y[i] - TIME - temp)
      # }
      # temp <- exp(tempest[1] + tempest[2] * data$Z1[count:N2] + tempest[3] * data$AZ1[count:N2] + tempest[4] * data$AZ2[count:N2])
      temp <- tempest[1] + tempest[2] * data$Z1[count:N2] + tempest[3] * data$AZ1[count:N2] + tempest[4] * data$AZ2[count:N2]
      alpha_mit[count:N2] <- data$delta[count:N2]*data$Sc[count:N2]*(data$Y[count:N2] - TIME_seq[o] - temp)
      
      # Second compute alpha.Exx and alpha.Exz
      tempxx <- 0; tempzx <- rep(0, 3)
      temp <- 0
      
      # for (i in count:N2){ # Y_i > t
      #   temp <- exp(tempest[1] + tempest[2] * data$Z1[i] + tempest[3] * data$AZ1[i] + tempest[4] * data$AZ2[i])
      #   
      #   tempxx <- tempxx + temp*data$delta[i]*data$Sc[i]
      #   
      #   tempzx <- tempzx + temp*data$delta[i]*data$Sc[i]*c(data$Z1[i], data$AZ1[i], data$AZ2[i])
      # }
      # temp <- exp(tempest[1] + tempest[2] * data$Z1[count:N2] + tempest[3] * data$AZ1[count:N2] + tempest[4] * data$AZ2[count:N2])
      temp <- tempest[1] + tempest[2] * data$Z1[count:N2] + tempest[3] * data$AZ1[count:N2] + tempest[4] * data$AZ2[count:N2]
      tempxx <- sum(temp*data$delta[count:N2]*data$Sc[count:N2])
      tempzx[1] <- sum(temp*data$delta[count:N2]*data$Sc[count:N2]*data$Z1[count:N2])
      tempzx[2] <- sum(temp*data$delta[count:N2]*data$Sc[count:N2]*data$AZ1[count:N2])
      tempzx[3] <- sum(temp*data$delta[count:N2]*data$Sc[count:N2]*data$AZ2[count:N2])
      
      alpha_Exx <- tempxx/N2
      
      alpha_Ezx <- tempzx/N2
      
      # Third compute alpha.q[j]
      tempExx <- alpha_Exx
      tempExxInv <- 1/tempExx
      for (j in 1:N2){
        tempxx <- 0
        # for (i in 1:N2){
        #   tempxx <- tempxx + alpha_mit[i]*phihat[j, i]*data$Sc[i]
        # }
        tempxx <- sum(alpha_mit[1:N2]*phihat[j, ]*data$Sc)
        alpha_q[j] <- tempxx/N2
        
        # Now everything is ready.
        
        tempxx <- alpha_mit[j] - alpha_q[j] - alpha_Ezx %*% t(xi[j, ] %*% solve(A))
        
        phi_alpha[j, o] <- tempExxInv*tempxx
      }
    }
    
  }
  
  return(list(phi_beta = phi_beta, phi_alpha = phi_alpha))
}

expit <- function(x) exp(x)/(1 + exp(x))

inverse <- function (f, lower = - 10, upper = 10) {
  function (y) {
    uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
  }
}

ESTSTD <- function(Y, X.PS, Xclassify, delta, eta_hat, treat, PS, Sc, theta_hat, t){
  N <- length(Y)
  xi <- rep(NA, N)
  zeta <- rep(NA, N)
  D1t <- matrix(NA, nrow = nrow(X.PS) + 1, ncol = N)
  treat_hat <- as.numeric(as.matrix(cbind(rep(1, N), Xclassify)) %*% eta_hat > 0)
  w <- (as.numeric(treat)*treat_hat + (1 - as.numeric(treat))*(1 - treat_hat))/PS
  
  
  pi_theta <- matrix(NA, nrow = N, ncol = ncol(X.PS))                         
                                                                                                                    
  #∂pi_i/∂theta                    
  for (i in 1:ncol(X.PS)){                                                                                                                              
    pi_theta[, i] <- exp(theta_hat %*% t(X.PS))*X.PS[, i]/(exp(theta_hat %*% t(X.PS)) + 1)^2
  }
  #∂w_i/∂pi_i
  treat_eta <- as.numeric(as.matrix(cbind(rep(1, N), Xclassify)) %*% eta_hat > 0) #决策所分配的治疗方案
  w_pi <- - (as.numeric(treat)*treat_eta + (1 - as.numeric(treat))*(1 - treat_eta))*(2*as.numeric(treat) - 1)/((2*as.numeric(treat) - 1)*PS + (1 - as.numeric(treat)))^2
  #D1
  D1 <- (apply(drop(w_pi)*pi_theta*drop(delta)*drop(Sc)*drop((Y > t))*drop(Y - t), 2, sum)*sum(w*delta*Sc*(Y > t)) - sum(w*delta*Sc*(Y > t)*(Y - t))*apply(drop(w_pi)*pi_theta*drop(delta)*drop(Sc)*drop((Y > t)), 2, sum))/sum(w*delta*Sc*(Y > t))^2
  #Phi1
  Phi1 <- PHI1(X.PS, theta_hat, treat)
  #D2
  D2 <- (sum(w*delta*Sc*(Y > t)*(Y - t))*(w*delta*(Y > t)*Sc^2) - (w*delta*Sc^2*(Y > t)*(Y - t))*sum(w*delta*Sc*(Y > t)))/sum(w*delta*Sc*(Y > t))^2
  #sum_phi3
  Phi3_result <- PHI3(Y, delta)
  Phi_3 <- Phi3_result$Phi3
  sum_phi3 <- Phi3_result$sum_phi3
  
  mIPWE_hat <- mean((Y > t)*w*delta*(Y - t)*Sc)
  SIPWE_hat <- mean((Y > t)*w*delta*Sc)
  MRL_IPWE <- mIPWE_hat/SIPWE_hat
  xi <- w*delta*Sc*((Y > t)*(Y - t) - (Y > t)*MRL_IPWE)/mean(w*delta*Sc*(Y > t))
  zeta <- D2 %*% t(Phi_3) + D1 %*% Phi1 + xi
  # zeta <- D2*sum_phi3 + D1 %*% Phi1 + xi
  
  sd <- sqrt(sum(zeta^2)/(N^2))
  # upper.conf <- MRL_IPWE + 1.96*sd
  # lower.conf <- MRL_IPWE - 1.96*sd
  
  # result <- list(zeta = zeta, sd = sd, upper.conf = upper.conf, lower.conf = lower.conf, xi = xi)
  result <- list(zeta = zeta, sd = sd, xi = xi)
  return(result)
}

ESTSTD_AIPW <- function(data, Y, X.PS, Xclassify, delta, eta_hat, treat, PS, Sc, theta_hat, t, alpha, tempest, TIME_seq, type){
  data <- data[order(Y), ]
  N <- length(Y)
  index <- which((Y[order(Y)] > t) == TRUE)[1]
  TT <- Y[order(Y)]; TT <- c(0, TT) 
  
  if (type == 'proportional'){
    treat_hat <- as.numeric(as.matrix(cbind(rep(1, N), Xclassify)) %*% eta_hat > 0)
    w <- (as.numeric(treat)*treat_hat + (1 - as.numeric(treat))*(1 - treat_hat))/PS
    
    m0 <- exp(cbind(data$X, data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% c(alpha[1], tempest[2:4]))
    mt <- exp(cbind(data$X, data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% c(alpha[index], tempest[2:4]))
    MRL_hat <- exp(cbind(data$X, data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% tempest)
    mu <- exp(data$X %*% t(alpha[1:index]) + cbind(data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% rbind(rep(tempest[2], index), rep(tempest[3], index), rep(tempest[4], index)))
    dT <- diff(TT[1:(index + 1)])
    temp <- (1/mu) %*% dT
    St <- exp(- temp)*m0/mt
    
    mA <- sum((Y > t)*w*delta*(Y - t)*Sc + (1 - w)*MRL_hat*St)
    SA <- sum((Y > t)*w*delta*Sc + (1 - w)*St)
    
    zeta_theta <- PHI1(X.PS, theta_hat, treat)
    
    zeta_Sc <- PHI3(Y, delta)
    
    phi_ab <- MRL_ESTSTD(data, alpha, TIME_seq, N, phihat = zeta_Sc$Phi3, type, tempest)
    zeta_ab <- phi_ab$phi_beta
    zeta_m0u <- phi_ab$phi_alpha
    
    #A
    A <- m0
    #B
    B <- mt
    #C
    C <- exp(- temp)
    
    #∂A/∂beta
    A_beta <- matrix(NA, nrow = N, ncol = 3)
    A_beta[, 1] <- exp(cbind(data$X, data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% c(alpha[1], tempest[2:4]))*data$Z1
    A_beta[, 2] <- exp(cbind(data$X, data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% c(alpha[1], tempest[2:4]))*treat_hat*data$Z1
    A_beta[, 3] <- exp(cbind(data$X, data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% c(alpha[1], tempest[2:4]))*treat_hat*data$Z2
    
    #∂B/∂beta
    B_beta <- matrix(NA, nrow = N, ncol = 3)
    B_beta[, 1] <- exp(cbind(data$X, data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% tempest)*data$Z1
    B_beta[, 2] <- exp(cbind(data$X, data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% tempest)*treat_hat*data$Z1
    B_beta[, 3] <- exp(cbind(data$X, data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% tempest)*treat_hat*data$Z2
    
    #∂C/∂beta
    C_beta <- matrix(NA, nrow = N, ncol = 3)
    mu_beta <- exp(data$X %*% t(alpha[1:index]) + 
                     cbind(data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% rbind(rep(tempest[2], index), rep(tempest[3], index), rep(tempest[4], index)))/
      (exp(data$X %*% t(alpha[1:index]) +
             cbind(data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% rbind(rep(tempest[2], index), rep(tempest[3], index), rep(tempest[4], index))))^2
    temp <- mu_beta %*% dT
    C_beta[, 1] <- C*temp*data$Z1
    C_beta[, 2] <- C*temp*treat_hat*data$Z1
    C_beta[, 3] <- C*temp*treat_hat*data$Z2
    
    #∂m/∂beta
    m_beta <- rep(NA, 3)
    m_beta[1] <- sum((1 - w)*(A_beta[, 1]*C + C_beta[, 1]*A))
    m_beta[2] <- sum((1 - w)*(A_beta[, 2]*C + C_beta[, 2]*A))
    m_beta[3] <- sum((1 - w)*(A_beta[, 3]*C + C_beta[, 3]*A))
    
    #∂S/∂beta
    S_beta <- rep(NA, 3)
    S_beta[1] <- sum((1 - w)*((A_beta[, 1]*C + C_beta[, 1]*A)*B - A*C*B_beta[, 1])/B^2)
    S_beta[2] <- sum((1 - w)*((A_beta[, 2]*C + C_beta[, 2]*A)*B - A*C*B_beta[, 2])/B^2)
    S_beta[3] <- sum((1 - w)*((A_beta[, 3]*C + C_beta[, 3]*A)*B - A*C*B_beta[, 3])/B^2)
    
    #D_ab
    D_ab <- rep(NA, 3)
    D_ab[1] <- (m_beta[1]*SA - S_beta[1]*mA)/SA^2
    D_ab[2] <- (m_beta[2]*SA - S_beta[2]*mA)/SA^2
    D_ab[3] <- (m_beta[3]*SA - S_beta[3]*mA)/SA^2
    
    # D_ab <- apply(D_ab, 2, sum)
    
    # D_m0u <- matrix(NA, nrow = N, ncol = index)
    D_m0u <- rep(NA, index)
    
    #∂MA/∂m0(t)
    #∂SA/∂m0(t)
    SA_m0t <- sum((1 - w)*exp(alpha[1])*C/exp(tempest[1])^2)
    D_m0u[index] <- - mA*SA_m0t/SA^2
    
    #∂MA/∂m0(Yi)
    mA_m0Yk <- rep(NA, index - 1)
    SA_m0Yk <- rep(NA, index - 1)
    for (k in 1:(index - 1)){
      #求∂mA/∂m0(Yk)
      mA_m0Yk[k] <- sum((1 - w)*m0*C*(TT[k + 1] - TT[k])/exp(alpha[k + 1])^2/exp(cbind(data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% c(tempest[2:4])))
      
      #∂SA/∂m0(Yk)
      SA_m0Yk[k] <- sum((1 - w)*exp(alpha[1])*C/exp(tempest[1])*(TT[k + 1] - TT[k])/exp(alpha[k + 1])^2/exp(cbind(data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% c(tempest[2:4])))
      D_m0u[k] <- (mA_m0Yk[k]*SA - mA*SA_m0Yk[k])/SA^2
    }
    
    
    # D_theta <- matrix(NA, nrow = N, ncol = ncol(X.PS))
    D_theta <- rep(NA, ncol(X.PS))
    #∂w/∂theta
    w_theta <- matrix(NA, nrow = N, ncol = ncol(X.PS))
    for (i in 1:ncol(X.PS)){
      w_theta[, i] <- (treat*treat_hat + (1 - treat)*treat_hat)*exp(theta_hat %*% t(X.PS))*X.PS[, i]*(exp(theta_hat %*% t(X.PS))*treat + (1 - treat))/(exp(theta_hat %*% t(X.PS))*treat + (1 - treat))^2 - (treat*treat_hat + (1 - treat)*treat_hat)*(1 + exp(theta_hat %*% t(X.PS)))*(treat*exp(theta_hat %*% t(X.PS))*X.PS[, i])/(exp(theta_hat %*% t(X.PS))*treat + (1 - treat))^2
    }
    
    #∂m/∂theta
    m_theta <- rep(NA, ncol(X.PS))
    for (i in 1:ncol(X.PS)){
      m_theta[i] <- sum(w_theta[, i]*delta*Sc*(Y > t)*(Y - t) - w_theta[, i]*B*St)
    }
    
    #∂S/∂theta
    S_theta <- rep(NA, ncol(X.PS))
    for (i in 1:ncol(X.PS)){
      S_theta[i] <- sum(w_theta[, i]*delta*Sc*(Y > t) - w_theta[, i]*St)
    }
    
    #D_theta
    for (i in 1:ncol(X.PS)){
      D_theta[i] <- (m_theta[i]*SA - S_theta[i]*mA)/SA^2
    }
    # D_theta <- apply(D_theta, 2, sum)
    
    D_Sc <- ((w*delta*(Y > t)*Sc^2)*mA - (w*delta*Sc^2*(Y > t)*(Y - t))*SA)/SA^2
    
    
    
    zeta1 <- D_theta %*% zeta_theta
    zeta23 <- D_ab %*% t(zeta_ab)[, - (N + 1)] 
    zeta4 <- D_m0u %*% t(phi_ab$phi_alpha[1:N, ]) 
    zeta5 <- D_Sc*zeta_Sc$sum_phi3
    zeta6 <- w*delta*Sc*(Y > t)*(Y - t - mA/SA)/(SA/N) + (1 - w)*(MRL_hat*St - St*mA/SA)/(SA/N)
    
    zeta <- as.numeric(zeta1) + as.numeric(zeta23) + as.numeric(zeta4) + as.numeric(zeta5) + as.numeric(zeta6)
    
    sd <- sqrt(sum(zeta^2)/N^2)
    
    result <- list(zeta = zeta, sd = sd, zeta1 = zeta1, zeta23 = zeta23, zeta4 = zeta4, zeta5 = zeta5, zeta6 = zeta6)
  }
  
  if (type == 'additive'){
    treat_hat <- as.numeric(as.matrix(cbind(rep(1, N), Xclassify)) %*% eta_hat > 0)
    w <- (as.numeric(treat)*treat_hat + (1 - as.numeric(treat))*(1 - treat_hat))/PS
    
    m0 <- cbind(data$X, data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% c(alpha[1], tempest[2:4])
    mt <- cbind(data$X, data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% c(alpha[index], tempest[2:4])
    MRL_hat <- cbind(data$X, data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% tempest
    mu <- data$X %*% t(alpha[1:index]) + cbind(data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% rbind(rep(tempest[2], index), rep(tempest[3], index), rep(tempest[4], index))
    dT <- diff(TT[1:(index + 1)])
    temp <- (1/mu) %*% dT
    St <- exp(- temp)*m0/mt
    
    mA <- sum((Y > t)*w*delta*(Y - t)*Sc + (1 - w)*MRL_hat*St)
    SA <- sum((Y > t)*w*delta*Sc + (1 - w)*St)
    
    zeta_theta <- PHI1(X.PS, theta_hat, treat)
    
    zeta_Sc <- PHI3(Y, delta)
    
    phi_ab <- MRL_ESTSTD(data, alpha, TIME_seq, N, phihat = zeta_Sc$Phi3, type, tempest)
    zeta_ab <- phi_ab$phi_beta
    zeta_m0u <- phi_ab$phi_alpha
    
    #A
    A <- m0
    #B
    B <- mt
    #C
    C <- exp(- temp)
    
    #∂A/∂beta
    A_beta <- matrix(NA, nrow = N, ncol = 3)
    A_beta[, 1] <- data$Z1
    A_beta[, 2] <- treat_hat*data$Z1
    A_beta[, 3] <- treat_hat*data$Z2
    
    #∂B/∂beta
    B_beta <- matrix(NA, nrow = N, ncol = 3)
    B_beta[, 1] <- data$Z1
    B_beta[, 2] <- treat_hat*data$Z1
    B_beta[, 3] <- treat_hat*data$Z2
    
    #∂C/∂beta
    C_beta <- matrix(NA, nrow = N, ncol = 3)
    mu_beta <- 1/(data$X %*% t(alpha[1:index]) + cbind(data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% rbind(rep(tempest[2], index), rep(tempest[3], index), rep(tempest[4], index)))^2
    temp <- mu_beta %*% dT
    C_beta[, 1] <- C*temp*data$Z1
    C_beta[, 2] <- C*temp*treat_hat*data$Z1
    C_beta[, 3] <- C*temp*treat_hat*data$Z2
    
    #∂m/∂beta
    m_beta <- rep(NA, 3)
    m_beta[1] <- sum((1 - w)*(A_beta[, 1]*C + C_beta[, 1]*A))
    m_beta[2] <- sum((1 - w)*(A_beta[, 2]*C + C_beta[, 2]*A))
    m_beta[3] <- sum((1 - w)*(A_beta[, 3]*C + C_beta[, 3]*A))
    
    #∂S/∂beta
    S_beta <- rep(NA, 3)
    S_beta[1] <- sum((1 - w)*((A_beta[, 1]*C + C_beta[, 1]*A)*B - A*C*B_beta[, 1])/B^2)
    S_beta[2] <- sum((1 - w)*((A_beta[, 2]*C + C_beta[, 2]*A)*B - A*C*B_beta[, 2])/B^2)
    S_beta[3] <- sum((1 - w)*((A_beta[, 3]*C + C_beta[, 3]*A)*B - A*C*B_beta[, 3])/B^2)
    
    #D_ab 
    D_ab <- rep(NA, 3)
    D_ab[1] <- (m_beta[1]*SA - S_beta[1]*mA)/SA^2
    D_ab[2] <- (m_beta[2]*SA - S_beta[2]*mA)/SA^2
    D_ab[3] <- (m_beta[3]*SA - S_beta[3]*mA)/SA^2
    
    # D_ab <- apply(D_ab, 2, sum)
    
    # D_m0u <- matrix(NA, nrow = N, ncol = index)
    D_m0u <- rep(NA, index)
    
    #∂MA/∂m0(t)
    #∂SA/∂m0(t)
    SA_m0t <- sum((1 - w)*m0*C/mt^2)
    D_m0u[index] <- - mA*SA_m0t/SA^2
    
    #∂MA/∂m0(Yi)
    mA_m0Yk <- rep(NA, index - 1)
    SA_m0Yk <- rep(NA, index - 1)
    for (k in 1:(index - 1)){
      #∂mA/∂m0(Yk)
      mA_m0Yk[k] <- sum((1 - w)*m0*C*(TT[k + 1] - TT[k])/(alpha[k + 1] + cbind(data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% c(tempest[2:4]))^2)
      
      #∂SA/∂m0(Yk)
      SA_m0Yk[k] <- sum((1 - w)*m0*C/mt*(TT[k + 1] - TT[k])/(alpha[k + 1] + cbind(data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% c(tempest[2:4]))^2)
      D_m0u[k] <- (mA_m0Yk[k]*SA - mA*SA_m0Yk[k])/SA^2
    }
    
    
    # D_theta <- matrix(NA, nrow = N, ncol = ncol(X.PS))
    D_theta <- rep(NA, ncol(X.PS))
    #∂w/∂theta
    w_theta <- matrix(NA, nrow = N, ncol = ncol(X.PS))
    for (i in 1:ncol(X.PS)){
      w_theta[, i] <- (treat*treat_hat + (1 - treat)*treat_hat)*exp(theta_hat %*% t(X.PS))*X.PS[, i]*(exp(theta_hat %*% t(X.PS))*treat + (1 - treat))/(exp(theta_hat %*% t(X.PS))*treat + (1 - treat))^2 - (treat*treat_hat + (1 - treat)*treat_hat)*(1 + exp(theta_hat %*% t(X.PS)))*(treat*exp(theta_hat %*% t(X.PS))*X.PS[, i])/(exp(theta_hat %*% t(X.PS))*treat + (1 - treat))^2
    }
    
    #∂m/∂theta
    m_theta <- rep(NA, ncol(X.PS))
    for (i in 1:ncol(X.PS)){
      m_theta[i] <- sum(w_theta[, i]*delta*Sc*(Y > t)*(Y - t) - w_theta[, i]*B*St)
    }
    
    #∂S/∂theta
    S_theta <- rep(NA, ncol(X.PS))
    for (i in 1:ncol(X.PS)){
      S_theta[i] <- sum(w_theta[, i]*delta*Sc*(Y > t) - w_theta[, i]*St)
    }
    
    #D_theta
    for (i in 1:ncol(X.PS)){
      D_theta[i] <- (m_theta[i]*SA - S_theta[i]*mA)/SA^2
    }
    # D_theta <- apply(D_theta, 2, sum)
    
    D_Sc <- ((w*delta*(Y > t)*Sc^2)*mA - (w*delta*Sc^2*(Y > t)*(Y - t))*SA)/SA^2
    
    
    
    zeta1 <- D_theta %*% zeta_theta
    zeta23 <- D_ab %*% t(zeta_ab)[, - 501]
    zeta4 <- D_m0u %*% t(phi_ab$phi_alpha[1:500, ])
    zeta5 <- D_Sc*zeta_Sc$sum_phi3
    zeta6 <- w*delta*Sc*(Y > t)*(Y - t - mA/SA)/(SA/N) + (1 - w)*(MRL_hat*St - St*mA/SA)/(SA/N)
    
    zeta <- as.numeric(zeta1) + as.numeric(zeta23) + as.numeric(zeta4) + as.numeric(zeta5) + as.numeric(zeta6)
    
    sd <- sqrt(sum(zeta^2)/N^2)
    
    result <- list(zeta = zeta, sd = sd, zeta1 = zeta1, zeta23 = zeta23, zeta4 = zeta4, zeta5 = zeta5, zeta6 = zeta6)
  }
  
  return(result)
}

Smooth_ESTSTD <- function(Y, X.PS, Xclassify, delta, eta_hat, treat, PS, Sc, theta_hat, t){
  N <- length(Y)
  xi <- rep(NA, N)
  zeta <- rep(NA, N)
  D1t <- matrix(NA, nrow = nrow(X.PS) + 1, ncol = N)
  
  c0 <- 4^(1/3)
  etaX <- as.matrix(cbind(rep(1, N), Xclassify)) %*% eta_hat
  h <- c0*N^(- 1/3)*sd(etaX)
  treat_hat <- pnorm(etaX/h, mean = 0, sd = 1)
  w <- (as.numeric(treat)*treat_hat + (1 - as.numeric(treat))*(1 - treat_hat))/PS
  
  
  pi_theta <- matrix(NA, nrow = N, ncol = ncol(X.PS))                         
  
  #∂pi_i/∂theta                    
  for (i in 1:ncol(X.PS)){                                                                                                                              
    pi_theta[, i] <- exp(theta_hat %*% t(X.PS))*X.PS[, i]/(exp(theta_hat %*% t(X.PS)) + 1)^2
  }
  #∂w_i/∂pi_i
  w_pi <- - (as.numeric(treat)*treat_hat + (1 - as.numeric(treat))*(1 - treat_hat))*(2*as.numeric(treat) - 1)/((2*as.numeric(treat) - 1)*PS + (1 - as.numeric(treat)))^2
  #D1_i
  D1 <- (apply(drop(w_pi)*pi_theta*drop(delta)*drop(Sc)*drop((Y > t))*drop(Y - t), 2, sum)*sum(w*delta*Sc*(Y > t)) - sum(w*delta*Sc*(Y > t)*(Y - t))*apply(drop(w_pi)*pi_theta*drop(delta)*drop(Sc)*drop((Y > t)), 2, sum))/sum(w*delta*Sc*(Y > t))^2
  #Phi1
  Phi1 <- PHI1(X.PS, theta_hat, treat)
  #D2_i
  D2 <- (sum(w*delta*Sc*(Y > t)*(Y - t))*(w*delta*(Y > t)*Sc^2) - (w*delta*Sc^2*(Y > t)*(Y - t))*sum(w*delta*Sc*(Y > t)))/sum(w*delta*Sc*(Y > t))^2
  #sum_phi3
  Phi3_result <- PHI3(Y, delta)
  Phi_3 <- Phi3_result$Phi3
  sum_phi3 <- Phi3_result$sum_phi3
  
  mIPWE_hat <- mean((Y > t)*w*delta*(Y - t)*Sc)
  SIPWE_hat <- mean((Y > t)*w*delta*Sc)
  MRL_IPWE <- mIPWE_hat/SIPWE_hat
  xi <- w*delta*Sc*((Y > t)*(Y - t) - (Y > t)*MRL_IPWE)/mean(w*delta*Sc*(Y > t))
  zeta <- as.numeric(D2)*sum_phi3 + D1 %*% Phi1 + as.numeric(xi)
  
  sd <- sqrt(sum(zeta^2)/(N^2))
  # upper.conf <- MRL_IPWE + 1.96*sd
  # lower.conf <- MRL_IPWE - 1.96*sd
  
  # result <- list(zeta = zeta, sd = sd, upper.conf = upper.conf, lower.conf = lower.conf, xi = xi)
  result <- list(zeta = zeta, sd = sd, xi = xi)
  return(result)
}

Smooth_ESTSTD_AIPW <- function(data, Y, X.PS, Xclassify, delta, eta_hat, treat, PS, Sc, theta_hat, t, alpha, tempest, TIME_seq, type){
  data <- data[order(Y), ]
  N <- length(Y)
  index <- which((Y[order(Y)] > t) == TRUE)[1]
  TT <- Y[order(Y)]; TT <- c(0, TT) 
  
  c0 <- 4^(1/3)
  etaX <- as.matrix(cbind(rep(1, N), Xclassify)) %*% eta_hat
  h <- c0*N^(- 1/3)*sd(etaX)
  treat_hat <- pnorm(etaX/h, mean = 0, sd = 1)
  w <- (as.numeric(treat)*treat_hat + (1 - as.numeric(treat))*(1 - treat_hat))/PS
  
  
  if (type == 'proportional'){
    
    m0 <- exp(cbind(data$X, data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% c(alpha[1], tempest[2:4]))
    mt <- exp(cbind(data$X, data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% c(alpha[index], tempest[2:4]))
    MRL_hat <- exp(cbind(data$X, data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% tempest)
    mu <- exp(data$X %*% t(alpha[1:index]) + cbind(data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% rbind(rep(tempest[2], index), rep(tempest[3], index), rep(tempest[4], index)))
    dT <- diff(TT[1:(index + 1)])
    temp <- (1/mu) %*% dT
    St <- exp(- temp)*m0/mt
    
    mA <- sum((Y > t)*w*delta*(Y - t)*Sc + (1 - w)*MRL_hat*St)
    SA <- sum((Y > t)*w*delta*Sc + (1 - w)*St)
    
    zeta_theta <- PHI1(X.PS, theta_hat, treat)
    
    zeta_Sc <- PHI3(Y, delta)
    
    phi_ab <- MRL_ESTSTD(data, alpha, TIME_seq, N, phihat = zeta_Sc$Phi3, type, tempest)
    zeta_ab <- phi_ab$phi_beta
    zeta_m0u <- phi_ab$phi_alpha
    
    #A
    A <- m0
    #B
    B <- mt
    #C
    C <- exp(- temp)
    
    #∂A/∂beta
    A_beta <- matrix(NA, nrow = N, ncol = 3)
    A_beta[, 1] <- exp(cbind(data$X, data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% c(alpha[1], tempest[2:4]))*data$Z1
    A_beta[, 2] <- exp(cbind(data$X, data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% c(alpha[1], tempest[2:4]))*treat_hat*data$Z1
    A_beta[, 3] <- exp(cbind(data$X, data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% c(alpha[1], tempest[2:4]))*treat_hat*data$Z2
    
    #∂B/∂beta
    B_beta <- matrix(NA, nrow = N, ncol = 3)
    B_beta[, 1] <- exp(cbind(data$X, data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% tempest)*data$Z1
    B_beta[, 2] <- exp(cbind(data$X, data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% tempest)*treat_hat*data$Z1
    B_beta[, 3] <- exp(cbind(data$X, data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% tempest)*treat_hat*data$Z2
    
    #∂C/∂beta
    C_beta <- matrix(NA, nrow = N, ncol = 3)
    mu_beta <- exp(data$X %*% t(alpha[1:index]) + 
                     cbind(data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% rbind(rep(tempest[2], index), rep(tempest[3], index), rep(tempest[4], index)))/
      (exp(data$X %*% t(alpha[1:index]) +
             cbind(data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% rbind(rep(tempest[2], index), rep(tempest[3], index), rep(tempest[4], index))))^2
    temp <- mu_beta %*% dT
    C_beta[, 1] <- C*temp*data$Z1
    C_beta[, 2] <- C*temp*treat_hat*data$Z1
    C_beta[, 3] <- C*temp*treat_hat*data$Z2
    
    #∂m/∂beta
    m_beta <- rep(NA, 3)
    m_beta[1] <- sum((1 - w)*(A_beta[, 1]*C + C_beta[, 1]*A))
    m_beta[2] <- sum((1 - w)*(A_beta[, 2]*C + C_beta[, 2]*A))
    m_beta[3] <- sum((1 - w)*(A_beta[, 3]*C + C_beta[, 3]*A))
    
    #∂S/∂beta
    S_beta <- rep(NA, 3)
    S_beta[1] <- sum((1 - w)*((A_beta[, 1]*C + C_beta[, 1]*A)*B - A*C*B_beta[, 1])/B^2)
    S_beta[2] <- sum((1 - w)*((A_beta[, 2]*C + C_beta[, 2]*A)*B - A*C*B_beta[, 2])/B^2)
    S_beta[3] <- sum((1 - w)*((A_beta[, 3]*C + C_beta[, 3]*A)*B - A*C*B_beta[, 3])/B^2)
    
    #D_ab 
    D_ab <- rep(NA, 3)
    D_ab[1] <- (m_beta[1]*SA - S_beta[1]*mA)/SA^2
    D_ab[2] <- (m_beta[2]*SA - S_beta[2]*mA)/SA^2
    D_ab[3] <- (m_beta[3]*SA - S_beta[3]*mA)/SA^2
    
    # D_ab <- apply(D_ab, 2, sum)
    
    # D_m0u <- matrix(NA, nrow = N, ncol = index)
    D_m0u <- rep(NA, index)
    
    #∂MA/∂m0(t)
    #∂SA/∂m0(t)
    SA_m0t <- sum((1 - w)*exp(alpha[1])*C/exp(tempest[1])^2)
    D_m0u[index] <- - mA*SA_m0t/SA^2
    
    #∂MA/∂m0(Yi)
    mA_m0Yk <- rep(NA, index - 1)
    SA_m0Yk <- rep(NA, index - 1)
    for (k in 1:(index - 1)){
      #∂mA/∂m0(Yk)
      mA_m0Yk[k] <- sum((1 - w)*m0*C*(TT[k + 1] - TT[k])/exp(alpha[k + 1])^2/exp(cbind(data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% c(tempest[2:4])))
      
      #∂SA/∂m0(Yk)
      SA_m0Yk[k] <- sum((1 - w)*exp(alpha[1])*C/exp(tempest[1])*(TT[k + 1] - TT[k])/exp(alpha[k + 1])^2/exp(cbind(data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% c(tempest[2:4])))
      D_m0u[k] <- (mA_m0Yk[k]*SA - mA*SA_m0Yk[k])/SA^2
    }
    
    
    # D_theta <- matrix(NA, nrow = N, ncol = ncol(X.PS))
    D_theta <- rep(NA, ncol(X.PS))
    #∂w/∂theta
    w_theta <- matrix(NA, nrow = N, ncol = ncol(X.PS))
    for (i in 1:ncol(X.PS)){
      w_theta[, i] <- as.numeric(treat*treat_hat + (1 - treat)*treat_hat)*exp(theta_hat %*% t(X.PS))*X.PS[, i]*(exp(theta_hat %*% t(X.PS))*treat + (1 - treat))/(exp(theta_hat %*% t(X.PS))*treat + (1 - treat))^2 - as.numeric(treat*treat_hat + (1 - treat)*treat_hat)*(1 + exp(theta_hat %*% t(X.PS)))*(treat*exp(theta_hat %*% t(X.PS))*X.PS[, i])/(exp(theta_hat %*% t(X.PS))*treat + (1 - treat))^2
    }
    
    #∂m/∂theta
    m_theta <- rep(NA, ncol(X.PS))
    for (i in 1:ncol(X.PS)){
      m_theta[i] <- sum(w_theta[, i]*delta*Sc*(Y > t)*(Y - t) - w_theta[, i]*B*St)
    }
    
    #∂S/∂theta
    S_theta <- rep(NA, ncol(X.PS))
    for (i in 1:ncol(X.PS)){
      S_theta[i] <- sum(w_theta[, i]*delta*Sc*(Y > t) - w_theta[, i]*St)
    }
    
    #D_theta
    for (i in 1:ncol(X.PS)){
      D_theta[i] <- (m_theta[i]*SA - S_theta[i]*mA)/SA^2
    }
    # D_theta <- apply(D_theta, 2, sum)
    
    D_Sc <- ((w*delta*(Y > t)*Sc^2)*mA - (w*delta*Sc^2*(Y > t)*(Y - t))*SA)/SA^2
    
    
    
    zeta1 <- D_theta %*% zeta_theta
    zeta23 <- D_ab %*% t(zeta_ab)[, - (N + 1)]
    zeta4 <- D_m0u %*% t(phi_ab$phi_alpha[1:N, ])
    zeta5 <- D_Sc*zeta_Sc$sum_phi3
    zeta6 <- w*delta*Sc*(Y > t)*(Y - t - mA/SA)/(SA/N) + (1 - w)*(MRL_hat*St - St*mA/SA)/(SA/N)
    
    zeta <- as.numeric(zeta1) + as.numeric(zeta23) + as.numeric(zeta4) + as.numeric(zeta5) + as.numeric(zeta6)
    
    sd <- sqrt(sum(zeta^2)/N^2)
    
    result <- list(zeta = zeta, sd = sd, zeta1 = zeta1, zeta23 = zeta23, zeta4 = zeta4, zeta5 = zeta5, zeta6 = zeta6)
  }
  
  if (type == 'additive'){
    
    m0 <- cbind(data$X, data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% c(alpha[1], tempest[2:4])
    mt <- cbind(data$X, data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% c(alpha[index], tempest[2:4])
    MRL_hat <- cbind(data$X, data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% tempest
    mu <- data$X %*% t(alpha[1:index]) + cbind(data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% rbind(rep(tempest[2], index), rep(tempest[3], index), rep(tempest[4], index))
    dT <- diff(TT[1:(index + 1)])
    temp <- (1/mu) %*% dT
    St <- exp(- temp)*m0/mt
    
    mA <- sum((Y > t)*w*delta*(Y - t)*Sc + (1 - w)*MRL_hat*St)
    SA <- sum((Y > t)*w*delta*Sc + (1 - w)*St)
    
    zeta_theta <- PHI1(X.PS, theta_hat, treat)
    
    zeta_Sc <- PHI3(Y, delta)
    
    phi_ab <- MRL_ESTSTD(data, alpha, TIME_seq, N, phihat = zeta_Sc$Phi3, type, tempest)
    zeta_ab <- phi_ab$phi_beta
    zeta_m0u <- phi_ab$phi_alpha
    
    #A
    A <- m0
    #B
    B <- mt
    #C
    C <- exp(- temp)
    
    #∂A/∂beta
    A_beta <- matrix(NA, nrow = N, ncol = 3)
    A_beta[, 1] <- data$Z1
    A_beta[, 2] <- treat_hat*data$Z1
    A_beta[, 3] <- treat_hat*data$Z2
    
    #∂B/∂beta
    B_beta <- matrix(NA, nrow = N, ncol = 3)
    B_beta[, 1] <- data$Z1
    B_beta[, 2] <- treat_hat*data$Z1
    B_beta[, 3] <- treat_hat*data$Z2
    
    #∂C/∂beta
    C_beta <- matrix(NA, nrow = N, ncol = 3)
    mu_beta <- 1/(data$X %*% t(alpha[1:index]) + cbind(data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% rbind(rep(tempest[2], index), rep(tempest[3], index), rep(tempest[4], index)))^2
    temp <- mu_beta %*% dT
    C_beta[, 1] <- C*temp*data$Z1
    C_beta[, 2] <- C*temp*treat_hat*data$Z1
    C_beta[, 3] <- C*temp*treat_hat*data$Z2
    
    #∂m/∂beta
    m_beta <- rep(NA, 3)
    m_beta[1] <- sum((1 - w)*(A_beta[, 1]*C + C_beta[, 1]*A))
    m_beta[2] <- sum((1 - w)*(A_beta[, 2]*C + C_beta[, 2]*A))
    m_beta[3] <- sum((1 - w)*(A_beta[, 3]*C + C_beta[, 3]*A))
    
    #∂S/∂beta
    S_beta <- rep(NA, 3)
    S_beta[1] <- sum((1 - w)*((A_beta[, 1]*C + C_beta[, 1]*A)*B - A*C*B_beta[, 1])/B^2)
    S_beta[2] <- sum((1 - w)*((A_beta[, 2]*C + C_beta[, 2]*A)*B - A*C*B_beta[, 2])/B^2)
    S_beta[3] <- sum((1 - w)*((A_beta[, 3]*C + C_beta[, 3]*A)*B - A*C*B_beta[, 3])/B^2)
    
    #D_ab
    D_ab <- rep(NA, 3)
    D_ab[1] <- (m_beta[1]*SA - S_beta[1]*mA)/SA^2
    D_ab[2] <- (m_beta[2]*SA - S_beta[2]*mA)/SA^2
    D_ab[3] <- (m_beta[3]*SA - S_beta[3]*mA)/SA^2
    
    # D_ab <- apply(D_ab, 2, sum)
    
    # D_m0u <- matrix(NA, nrow = N, ncol = index)
    D_m0u <- rep(NA, index)
    
    #∂MA/∂m0(t)
    #∂SA/∂m0(t)
    SA_m0t <- sum((1 - w)*m0*C/mt^2)
    D_m0u[index] <- - mA*SA_m0t/SA^2
    
    #∂MA/∂m0(Yi)
    mA_m0Yk <- rep(NA, index - 1)
    SA_m0Yk <- rep(NA, index - 1)
    for (k in 1:(index - 1)){
      #∂mA/∂m0(Yk)
      mA_m0Yk[k] <- sum((1 - w)*m0*C*(TT[k + 1] - TT[k])/(alpha[k + 1] + cbind(data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% c(tempest[2:4]))^2)
      
      #∂SA/∂m0(Yk)
      SA_m0Yk[k] <- sum((1 - w)*m0*C/mt*(TT[k + 1] - TT[k])/(alpha[k + 1] + cbind(data$Z1, treat_hat*data$Z1, treat_hat*data$Z2) %*% c(tempest[2:4]))^2)
      D_m0u[k] <- (mA_m0Yk[k]*SA - mA*SA_m0Yk[k])/SA^2
    }
    
    
    # D_theta <- matrix(NA, nrow = N, ncol = ncol(X.PS))
    D_theta <- rep(NA, ncol(X.PS))
    #∂w/∂theta
    w_theta <- matrix(NA, nrow = N, ncol = ncol(X.PS))
    for (i in 1:ncol(X.PS)){
      w_theta[, i] <- (treat*treat_hat + (1 - treat)*treat_hat)*exp(theta_hat %*% t(X.PS))*X.PS[, i]*(exp(theta_hat %*% t(X.PS))*treat + (1 - treat))/(exp(theta_hat %*% t(X.PS))*treat + (1 - treat))^2 - (treat*treat_hat + (1 - treat)*treat_hat)*(1 + exp(theta_hat %*% t(X.PS)))*(treat*exp(theta_hat %*% t(X.PS))*X.PS[, i])/(exp(theta_hat %*% t(X.PS))*treat + (1 - treat))^2
    }
    
    #∂m/∂theta
    m_theta <- rep(NA, ncol(X.PS))
    for (i in 1:ncol(X.PS)){
      m_theta[i] <- sum(w_theta[, i]*delta*Sc*(Y > t)*(Y - t) - w_theta[, i]*B*St)
    }
    
    #∂S/∂theta
    S_theta <- rep(NA, ncol(X.PS))
    for (i in 1:ncol(X.PS)){
      S_theta[i] <- sum(w_theta[, i]*delta*Sc*(Y > t) - w_theta[, i]*St)
    }
    
    #D_theta
    for (i in 1:ncol(X.PS)){
      D_theta[i] <- (m_theta[i]*SA - S_theta[i]*mA)/SA^2
    }
    # D_theta <- apply(D_theta, 2, sum)
    
    D_Sc <- ((w*delta*(Y > t)*Sc^2)*mA - (w*delta*Sc^2*(Y > t)*(Y - t))*SA)/SA^2
    
    
    
    zeta1 <- D_theta %*% zeta_theta
    zeta23 <- D_ab %*% t(zeta_ab)[, - 501]
    zeta4 <- D_m0u %*% t(phi_ab$phi_alpha[1:500, ])
    zeta5 <- D_Sc*zeta_Sc$sum_phi3
    zeta6 <- w*delta*Sc*(Y > t)*(Y - t - mA/SA)/(SA/N) + (1 - w)*(MRL_hat*St - St*mA/SA)/(SA/N)
    
    zeta <- as.numeric(zeta1) + as.numeric(zeta23) + as.numeric(zeta4) + as.numeric(zeta5) + as.numeric(zeta6)
    
    sd <- sqrt(sum(zeta^2)/N^2)
    
    result <- list(zeta = zeta, sd = sd, zeta1 = zeta1, zeta23 = zeta23, zeta4 = zeta4, zeta5 = zeta5, zeta6 = zeta6)
  }
  
  return(result)
}

PHI1 <- function(X, betaHat, treat){
  p <- length(betaHat)
  n <- nrow(X)
  y <- treat
  X <- as.matrix(X)
  PhiHat <- matrix(NA, nrow = p, ncol = n)
  UHat <- rep(NA, p)
  uHat <- matrix(NA, nrow = p, ncol = n)
  for (j in 1:p){
    for (m in 1:n){
      uHat[j, m] <- (y[m] - exp(X[m, ] %*% betaHat)/(1 + exp(X[m, ] %*% betaHat)))*X[m, j]
    }
  }
  Sigma <- matrix(NA, nrow = p, ncol = p)
  for (i in 1:p){
    for (j in 1:p){
      temp <- 0
      for (m in 1:n){
        temp <- temp - X[m, i]*X[m, j]*exp(X[m, ] %*% betaHat)/(1 + exp(X[m, ] %*% betaHat))^2
        Sigma[i, j] <- temp
        
      }
    }
  }
  Sigma <- solve(- Sigma/n)
  for (i in 1:n){
    PhiHat[, i] <- Sigma %*% uHat[, i]
  }
  return(PhiHat)
}

PHI3 <- function(Y, delta){
  if (sum(delta) == length(Y)){
    Phi <- matrix(NA, nrow = N, ncol = N)
    sum_phi3 <- apply(Phi, 2, mean)
  }else{
    KM <- survfit(Surv(Y, delta) ~ 1)
    N <- length(Y)
    S <- rep(NA, N)
    for (i in 1:N){
      temp <- which(Y[i] <= KM$time)[1]
      if (KM$surv[temp] == 0){
        temp <- temp - 1
      }
      S[i] <- KM$surv[temp] 
    }
    
    status <- 1 - delta
    KM_censor <- survfit(Surv(Y, status) ~ 1)
    Sc <- rep(NA, N)
    G <- rep(NA, N)
    for (i in 1:N){
      temp <- which(Y[i] <= KM_censor$time)[1]
      if (KM_censor$surv[temp] == 0){
        temp <- temp - 1
      }
      G[i] <-  KM_censor$surv[temp]
      Sc[i] <- 1/KM_censor$surv[temp]
    }
    
    Censor_time <- Y[delta == 0]
    Censor_time <- Censor_time[order(Censor_time)]
    Censor_position <- which(delta[order(Y)]==0)
    TT <- Y[order(Y)]; TT <- c(0, TT) 
    
    
    LAMBDA <- KM_censor$cumhaz; LAMBDA <- c(0, LAMBDA) 
    #Pi^c(t)=P(Yi>=t)
    pict <- S[order(Y)]
    #G(t)
    G <- G[order(Y)]
    
    N <- length(Y)
    
    ############Compute PHI###############
    #PHI
    Phi <- matrix(NA, nrow = N, ncol = N) 
    
    for (o in 1:length(Censor_position)){
      print(o)
      tau <- Censor_position[o]
      
      index <- which((Y[order(Y)] >= Censor_time[o]) == TRUE)[1]
      dLAMBDA <- matrix(rep(diff(LAMBDA[1:(index + 1)]), index), nrow = index, ncol = index)
      dLAMBDA[lower.tri(dLAMBDA)] <- 0
      
      
      YY <- matrix(rep(Y, index), nrow = N, ncol = index)
      TTT <- matrix(rep(TT[2:(index + 1)], N), nrow = N, ncol = index, byrow = T)
      jict <- as.matrix(YY >= TTT)
      
      
      mat_delta <- matrix(rep((1 - delta), index), nrow = N, ncol = index)
      temp <- as.matrix(YY <= TTT)
      nict <- mat_delta*temp
      rm(temp); rm(mat_delta); rm(TTT)
      
      Mict <- nict - jict %*% dLAMBDA
      rm(nict); rm(jict)
      Mict <- cbind(rep(0, N), Mict)
      dMict <- diff(t(Mict))
      rm(Mict); rm(dLAMBDA)
      
      Phi[, tau] <- - (1/pict[1:index]) %*% dMict
      Phi[, tau] <- Phi[, tau]*G[index]
      rm(dMict)
      
      
      if (o == 1){
        for (m in 1:(Censor_position[o] - 1)){
          Phi[, m] <- rep(0, N)
        }
      }
      
      if (o == length(Censor_position)){
        for (m in min(Censor_position[o] + 1, N):N){
          Phi[, m] <- Phi[, Censor_position[o]]
        }
      }else{
        for (m in (Censor_position[o] + 1):(Censor_position[o + 1] - 1)){
          Phi[, m] <- Phi[, Censor_position[o]]
        }
      }
      sum_phi3 <- apply(Phi, 1, mean)
    }
  }
  
  
  return(list(Phi3 = Phi, sum_phi3 = sum_phi3))
}

bootstrap_sd <- function(Y, Xclassify, delta, treat, w, Sc, PS, eta_hat){
  boot_V <- rep(0, length(Y))
  boot_SV <- rep(0, length(Y))
  for (i in 1:length(Y)){
    boot_index <- sample(1:length(Y), size = length(Y), replace = TRUE, prob = NULL)
    treat_eta <- as.numeric(as.matrix(cbind(rep(1, N), Xclassify[boot_index, ])) %*% eta_hat > 0) 
    w <- (as.numeric(treat[boot_index])*treat_eta + (1 - as.numeric(treat[boot_index]))*(1 - treat_eta))/PS[boot_index]
    temp1 <- (Y[boot_index] > t)*w[boot_index]*delta[boot_index]*(Y[boot_index] - t)*Sc[boot_index]
    Temp1 <- (Y[boot_index] > t)*w[boot_index]*delta[boot_index]*Sc[boot_index]
    boot_V[i] <- mean(temp1)/mean(Temp1)
    
    
    c0 <- 4^(1/3)
    etaX <- as.matrix(cbind(rep(1, N), Xclassify[boot_index, ])) %*% eta_hat
    h <- c0*N^(- 1/3)*sd(etaX)
    Streat_eta <- pnorm(etaX/h, mean = 0, sd = 1)
    Sw <- (as.numeric(treat[boot_index])*Streat_eta + (1 - as.numeric(treat[boot_index]))*(1 - Streat_eta))/PS[boot_index]
    temp2 <- (Y[boot_index] > t)*Sw[boot_index]*delta[boot_index]*(Y[boot_index] - t)*Sc[boot_index]
    Temp2 <- (Y[boot_index] > t)*Sw[boot_index]*delta[boot_index]*Sc[boot_index]
    boot_SV[i] <- mean(temp2)/mean(Temp2)
  }
  boot_std1 <- sd(boot_V)
  boot_std2 <- sd(boot_SV)
  return(list(boot_std1 = boot_std1, boot_std2 = boot_std2, boot_V = boot_V, boot_SV = boot_SV))
}

findrep <- function(a){
  q <- 0
  for (i in 1:(length(a)-1)){
    for (j in (i + 1):length(a)){
      if (a[i] == a[j]){
        q <- q + 1
      }
    }
  }
  return(q)
}