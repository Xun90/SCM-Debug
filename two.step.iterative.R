two.step.iterative <- function(Y1pre,Y0pre,X1.scaled,X0.scaled,SV){
  #SV - predictor weights defined by the user
  ##Solve non-Archimedean problem (8)
  Tpre <- dim(Y0pre)[1]
  nDonors <- dim(Y0pre)[2]
  #QP setup
  A <- matrix(rep(1,nDonors), ncol = nDonors)
  b <- 1
  r <- 0
  l <- matrix(rep(0,nDonors), nrow = nDonors)
  u <- matrix(rep(1,nDonors), nrow = nDonors)
  #Loop on 10 epsilon values (0.1^1 ... 0.1^10) to find the best performer 
  L_upper = matrix(0, 10, 2) #for storing the optimal value of the upper level problem
  L_lower = matrix(0, 10, 2) #for storing the optimal value of the lower level problem
  W_ipop = matrix(0, nDonors, 10) #for storing the optimal W weights from QP_Solver1
  W_lowr = matrix(0, nDonors, 10) #for storing the optimal W weights from QP_Solver2
  for (i in 1:10){
    eps <- 0.1^(i) #epsilon - penalty term
    c <- (-t(X0.scaled) %*% diag(SV1) %*% X1.scaled) - eps * t(Y0pre) %*% Y1pre
    H <- t(X0.scaled) %*% diag(SV1) %*% X0.scaled + eps * t(Y0pre) %*% Y0pre
    #run QP
    QP_ipop <- ipop(c = c, H = H, A = A, b = b, l = l, u = u, r = r, 
                    margin = 0.0005, maxiter = 1000, sigf = 7, bound = 10) #QP_Solver1
    QP_lowr <- LowRankQP(Vmat = H, dvec = c, Amat = A, bvec = b, uvec = u, 
                         method = "LU") #QP_Solver2
    W_ipop[,i] <- matrix(QP_ipop@primal, nrow = nDonors)
    W_lowr[,i] <- QP_lowr$alpha
    L_upper[i,1] <- 1/Tpre * t(Y1pre - Y0pre %*% W_ipop[,i]) %*% 
      (Y1pre - Y0pre %*% W_ipop[,i])
    L_upper[i,2] <- 1/Tpre * t(Y1pre - Y0pre %*% W_lowr[,i]) %*% 
      (Y1pre - Y0pre %*% W_lowr[,i])
    L_lower[i,1] <- t(X1.scaled - X0.scaled %*% W_ipop[,i]) %*% diag(SV1) %*% 
      (X1.scaled - X0.scaled %*% W_ipop[,i])
    L_lower[i,2] <- t(X1.scaled - X0.scaled %*% W_lowr[,i]) %*% diag(SV1) %*% 
      (X1.scaled - X0.scaled %*% W_lowr[,i])
  }
  ##Use the first step of the two-step procedure in Section 4.1 to determine epsilon
  c1 <- (-t(X0.scaled) %*% diag(SV) %*% X1.scaled)
  H1 <- t(X0.scaled) %*% diag(SV) %*% X0.scaled
  #run QP
  QP_ipop1 <- ipop(c = c1, H = H1, A = A, b = b, l = l, u = u, r = r, 
                   margin = 0.0005, maxiter = 1000, sigf = 7, bound = 10) #QP_Solver1
  QP_lowr1 <- LowRankQP(Vmat = H1, dvec = c1, Amat = A, bvec = b, uvec = u, 
                        method = "LU") #QP_Solver2
  W_ipop1 <- matrix(QP_ipop1@primal, nrow = nDonors)
  W_lowr1 <- QP_lowr1$alpha
  W1 <- cbind(W_ipop1, W_lowr1)
  obj_left <- t(X1.scaled) %*% diag(SV) %*% X1.scaled 
  Lw_ipop1 <- obj_left + 2 * (t(c1) %*% W_ipop1 + 1/2 * t(W_ipop1) %*% H1 %*% W_ipop1)
  Lw_lowr1 <- obj_left + 2 * (t(c1) %*% W_lowr1 + 1/2 * t(W_lowr1) %*% H1 %*% W_lowr1)
  Lw1 <- c(Lw_ipop1, Lw_lowr1)
  
  two.step.iterative.out <- list(W = cbind(W_ipop,W_lowr), W1 = W1, V = SV, 
                         L_upper = L_upper, L_lower = L_lower, Lw1 = Lw1)
  return(two.step.iterative.out)
}