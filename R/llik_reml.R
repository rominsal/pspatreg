## llik REML concentrated in rho-delta-phi.
llikc <- function(param, env) {
  if (env$type %in% c("sar", "sdm")) {
    rho <- param[1]
    delta <- NULL
    if (env$cor == "ar1") 
      phi <- param[2]
    else phi <- NULL
  } else if (env$type %in% c("sem", "sdem")) {
    delta <- param[1]
    rho <- NULL
    if (env$cor == "ar1") 
      phi <- param[2]
    else phi <- NULL
  } else if (env$type == "sarar") {
    rho <- param[1]
    delta <- param[2]
    if (env$cor == "ar1")
      phi <- param[3]
    else phi <- NULL
  } else {
    rho <- delta <- NULL
    if (env$cor == "ar1") 
      phi <- param[1]
    else phi <- NULL
  }
  nfull <- env$nfull
  nsp <- env$nsp
  nt <- env$nt
  if ((is.null(nfull) || is.null(nsp)) || is.null(nt))
    stop("nfull or nsp or nt is NULL")
  y <- env$y
  nfull <- length(y)
  X <- env$Xfull
  if (!is.null(env$Zfull)) {
    Z <- env$Zfull
  } else Z <- matrix(0, nrow = length(y), ncol = 1)
  D <- env$D
  sig2u <- env$sig2u
  Wsp <- env$Wsp
  Ifull <- Diagonal(nfull)
  Isp <- Diagonal(nsp)
  # Change It when there is temporal correlation...
  It <- Diagonal(nt)
  if (!is.null(rho)) {
    A1 <- Isp - rho*Wsp
    ldetA1 <- do_ldet(rho, env)
  } else {
    A1 <- Isp
    ldetA1 <- 0
  }
  if (!is.null(delta)) {
    A2 <- Isp - delta*Wsp
    ldetA2 <- do_ldet(delta, env)
  } else {
    A2 <- Isp
    ldetA2 <- 0
  }
  if (nt > 1) {
    if (!is.null(phi)) { 
      call_Omega <- build_Omega_ar1(phi, nt)
      Omega <- call_Omega$Omega
      Omegainv <- call_Omega$Omegainv
      LcholOmegainv <- chol(Omegainv)
    } else { # phi <- NULL
      Omega <- Omegainv <- LcholOmegainv <- It
    }
    # fast index nt, slow index nsp
    ystar <- matrix(RH(A2 %*% A1, 
                       RH(LcholOmegainv, 
                          array(y, dim = c(ncol(LcholOmegainv), 
                                           ncol(A1))))), 
                    ncol = 1)
    Xstar <- kronecker(A2, LcholOmegainv) %*% X
    Zstar <- kronecker(A2, LcholOmegainv) %*% Z
  } else { # nt == 1
    Omega <- Omegainv <- LcholOmegainv <- It
    ystar <- A2 %*% (A1 %*% y)
    Xstar <- A2 %*% X
    Zstar <- A2 %*% Z    
  }
  mat <- construct_matrices(Xstar, Zstar, ystar)
  if (env$nvarspt > 0 || env$nvarnopar > 0) {
    # MATRIX C IN (12). PAPER SAP
    C <- Matrix( construct_block(
      mat$XtX, 
      t(mat$ZtX*env$G_eff), 
      mat$ZtX, 
      t(mat$ZtZ*env$G_eff)))
    H <- (1/sig2u)*C + D
    Hinv <- solve(H)
    G <- Diagonal(length(env$G_eff), x = env$G_eff)
    ## Formulae to compute P, ldetV and ldetV.plus.ldetXtVinvX
    ## from paper Harville (1977) pp. 326
    P <- (1/sig2u)*Ifull - 
      (1/sig2u^2)*cbind(Xstar, Zstar %*% G) %*% Hinv %*% 
      t(cbind(Xstar, Zstar))
    ldetV <- nfull*log(sig2u) + 
      determinant(Diagonal(length(env$G_eff)) + 
                            (1/sig2u)*mat$ZtZ*env$G_eff)$modulus    
  }  else { # Only fixed effects
    C <- Matrix(mat$XtX)
    H <- (1/sig2u)*C 
    Hinv <- solve(H)
    P <- (1/sig2u)*Ifull - 
      (1/sig2u^2)*(Xstar %*% Hinv) %*% t(Xstar)
    ldetV <- nfull*log(sig2u)
  }
  ystar_P_ystar <- t(ystar) %*% (P %*% ystar)
  log_likc <- -0.5*(ldetV + ystar_P_ystar) + 
    nt*ldetA1 + nt*ldetA2
  return(as.numeric(-log_likc))
}
##########################################################
## llik REML concentrated in rho-delta-phi.
llikc_reml <- function(param, env) {
  if (env$type %in% c("sar", "sdm")) {
    rho <- param[1]
    delta <- NULL
    if (env$cor == "ar1") 
      phi <- param[2]
    else phi <- NULL
  } else if (env$type %in% c("sem", "sdem")) {
    delta <- param[1]
    rho <- NULL
    if (env$cor == "ar1") 
      phi <- param[2]
    else phi <- NULL
  } else if (env$type == "sarar") {
    rho <- param[1]
    delta <- param[2]
    if (env$cor == "ar1") 
      phi <- param[3]
    else phi <- NULL
  } else {
    rho <- delta <- NULL
    if (env$cor == "ar1") {
      phi <- param[1]
    } else phi <- NULL
  }
  nfull <- env$nfull
  nsp <- env$nsp
  nt <- env$nt
  if ((is.null(nfull) || is.null(nsp)) || is.null(nt))
    stop("nfull or nsp or nt is NULL")
  y <- env$y
  X <- env$Xfull
  if (!is.null(env$Zfull)) {
    Z <- env$Zfull
  } else Z <- matrix(0, nrow = length(y), ncol = 1)
  D <- env$D
  Wsp <- env$Wsp
  sig2u <- env$sig2u
  Ifull <- Diagonal(nfull)
  Isp <- Diagonal(nsp)
  # Change It when there is temporal correlation...
  It <- Diagonal(nt)
  if (!is.null(rho)) {
    A1 <- Isp - rho*Wsp
    ldetA1 <- do_ldet(rho, env)
  } else {
    A1 <- Isp
    ldetA1 <- 0
  }
  if (!is.null(delta)) {
    A2 <- Isp - delta*Wsp
    ldetA2 <- do_ldet(delta, env)
  } else {
    A2 <- Isp
    ldetA2 <- 0
  }
  if (nt > 1) {
    if (!is.null(phi)) { 
      call_Omega <- build_Omega_ar1(phi, nt)
      Omega <- call_Omega$Omega
      Omegainv <- call_Omega$Omegainv
      LcholOmegainv <- chol(Omegainv)
    } else { # phi <- NULL
      Omega <- Omegainv <- LcholOmegainv <- It
    }
    # fast index nt, slow index nsp
    ystar <- matrix(RH(A2 %*% A1, 
                       RH(LcholOmegainv, 
                          array(y, dim = c(ncol(LcholOmegainv), 
                                           ncol(A1))))), 
                    ncol = 1)
    Xstar <- kronecker(A2, LcholOmegainv) %*% X
    Zstar <- kronecker(A2, LcholOmegainv) %*% Z
  } else { # nt == 1
    Omega <- Omegainv <- LcholOmegainv <- It
    ystar <- A2 %*% (A1 %*% y)
    Xstar <- A2 %*% X
    Zstar <- A2 %*% Z    
  }
  mat <- construct_matrices(Xstar, Zstar, ystar) 
  if (env$nvarspt > 0 || env$nvarnopar > 0) {
    # MATRIX C IN (12). PAPER SAP
    C <- Matrix( 
      construct_block(mat$XtX, 
                      t(mat$ZtX*env$G_eff), 
                      mat$ZtX, 
                      t(mat$ZtZ*env$G_eff)))
    H <- (1/sig2u)*C + D
    Hinv <- try(solve(H))
    if (inherits(Hinv, "try-error")) 
      Hinv <- ginv(as.matrix(H))
    G <- Diagonal(length(env$G_eff), x = env$G_eff)
    ## Formulae to compute P, ldetV and ldetV.plus.ldetXtVinvX
    ## from paper Harville (1977) pp. 326
    P <- (1/sig2u)*Ifull - 
      (1/sig2u^2)*cbind(Xstar, Zstar %*% G) %*% Hinv %*% 
      t(cbind(Xstar, Zstar))
  } else { # Only fixed effects
    C <- Matrix(mat$XtX)
    H <- (1/sig2u)*C 
    Hinv <- solve(H)
    browser()
    P <- (1/sig2u)*Ifull - 
      (1/sig2u^2)*(Xstar %*% Hinv) %*% t(Xstar)
  }
  ldetV.plus.ldetXtVinvX <- nfull*log(sig2u) + 
    determinant(H)$modulus
  ystar_P_ystar <- t(ystar) %*% (P %*% ystar)
  log_likc_reml <- -0.5*(ldetV.plus.ldetXtVinvX + 
                           ystar_P_ystar) + 
                    nt*ldetA1 + nt*ldetA2
  return(as.numeric(-log_likc_reml))
}
## analytic score llik REML 2d concentrated in rho.
## REPASAR CAMBIAR PARA 2D Y 3D Y AÑADIR CORRELACIÓN... 
ansco_llikc_reml_2d <- function(param, env) {
  if (env$type == "sar") {
    rho <- param[1]
    delta <- NULL
    if (env$cor == "ar1") {
      phi <- param[2]
    } else phi <- NULL
  } else if (env$type %in% c("sem", "sdem")) {
    delta <- param[1]
    rho <- NULL
    if (env$cor == "ar1") 
      phi <- param[2]
    else phi <- NULL
  } else if (env$type == "sarar") {
    rho <- param[1]
    delta <- param[2]
    if (env$cor == "ar1") {
      phi <- param[3]
    } else phi <- NULL
  } else {
    rho <- delta <- NULL
    if (env$cor == "ar1") {
      phi <- param[1]
    } else phi <- NULL
  }
  y <- env$y
  nfull <- length(y)
  X <- env$Xfull
  if (!is.null(env$Zfull)) {
    Z <- env$Zfull
  } else Z <- matrix(0, nrow = length(y), ncol = 1)
  D <- env$D
  Wsp <- env$Wsp
  if (!is.null(Wsp)) nsp <- nrow(Wsp) else nsp <- NULL
  sig2u <- env$sig2u
  In <- Diagonal(nfull)
  if (!is.null(rho)) {
    A1 <- In - rho*Wsp } else { A1 <- In }
  if (!is.null(delta)) { 
    A2 <- In - delta*Wsp } else  { A2 <- In }
  ystar <- Matrix(A2 %*% (A1 %*% y))
  Xstar <- Matrix(A2 %*% X)
  Zstar <- Matrix(A2 %*% Z)
  ## CALCULAR ldetV, P, ldetV.plus.ldetXtVinvX
  mat <- construct_matrices(Xstar, Zstar, ystar)
  if (env$nvarspt > 0 || env$nvarnopar > 0) {
    # MATRIX C IN (12). PAPER SAP
    C <- Matrix( 
      construct_block(mat$XtX, 
                      t(mat$ZtX*env$G_eff), 
                      mat$ZtX, 
                      t(mat$ZtZ*env$G_eff)))
    H <- (1/sig2u)*C + D
    Hinv <- solve(H)
    G <- Diagonal(length(env$G_eff), x = env$G_eff)
    ## Formulae to compute P, ldetV and ldetV.plus.ldetXtVinvX
    ## from paper Harville (1977) pp. 326
    P <- (1/sig2u)*In - 
      (1/sig2u^2)*cbind(Xstar, Zstar %*% G) %*% Hinv %*% 
      t(cbind(Xstar, Zstar))
  } else {  # Only fixed effects
    C <- Matrix(mat$XtX)
    H <- (1/sig2u)*C 
    Hinv <- solve(H)
    P <- (1/sig2u)*In - 
      (1/sig2u^2)*(Xstar %*% Hinv) %*% t(Xstar)
  }
  if (env$type == "sar") {
    Wsp_y <- Matrix(Wsp %*% y)
    A1invWsp <- solve(A1, Wsp)
    score_REML <- as.numeric(t(P %*% ystar) %*% Wsp_y -
      sum(diag(as.matrix(A1invWsp))))
  } else { score_REML <- NULL }# FALTA IMPLEMENTAR SEM Y SARAR
  return(score_REML)
}  
## analytic hessian llik REML 2d concentrated in rho.
## REPASAR. CAMBIAR PARA 2D Y 3D Y AÑADIR CORRELACIÓN
anhess_llikc_reml_2d <- function(param, env) {
  if (env$type == "sar") {
    rho <- param[1]
    delta <- NULL
    if (env$cor == "ar1") {
      phi <- param[2]
    } else phi <- NULL
  } else if (env$type == "sem") {
    delta <- param[1]
    rho <- NULL
    if (env$cor == "ar1") {
      phi <- param[2]
    } else phi <- NULL
  } else if (env$type == "sarar") {
    rho <- param[1]
    delta <- param[2]
    if (env$cor == "ar1") {
      phi <- param[3]
    } else phi <- NULL
  } else {
    rho <- delta <- NULL
    if (env$cor == "ar1") {
      phi <- param[1]
    } else phi <- NULL
  }
  y <- env$y
  X <- env$Xfull
  if (!is.null(env$Zfull)) {
    Z <- env$Zfull
  } else Z <- matrix(0, nrow = length(y), ncol = 1)
  D <- env$D
  Wsp <- env$Wsp
  nsp <- nrow(Wsp)
  sig2u <- env$sig2u
  In <- Diagonal(nsp)
  if (!is.null(rho)) { 
    A1 <- In - rho*Wsp } else { A1 <- In }
  if (!is.null(delta)) { 
    A2 <- In - delta*Wsp } else { A2 <- In }
  ystar <- Matrix(A2 %*% (A1 %*% y))
  Xstar <- Matrix(A2 %*% X)
  Zstar <- Matrix(A2 %*% Z)
  ## CALCULAR ldetV, P, ldetV.plus.ldetXtVinvX
  mat <- construct_matrices(Xstar, Zstar, ystar)
  if (env$nvarspt > 0 || env$nvarnopar > 0) {
    C <- Matrix( 
      construct_block(mat$XtX, 
                      t(mat$ZtX*env$G_eff), 
                      mat$ZtX, 
                      t(mat$ZtZ*env$G_eff)))
    H <- (1/sig2u)*C + D
    Hinv <- solve(H)
    G <- Diagonal(length(env$G_eff), x = env$G_eff)
    ## Formulae to compute P, ldetV and ldetV.plus.ldetXtVinvX
    ## from paper Harville (1977) pp. 326
    P <- (1/sig2u)*In - 
      (1/sig2u^2)*cbind(Xstar, Zstar %*% G) %*% Hinv %*% 
      t(cbind(Xstar, Zstar))
  } else { # Only fixed effects
    C <- Matrix(mat$XtX)
    H <- (1/sig2u)*C
    Hinv <- solve(H)
    P <- (1/sig2u)*In - 
      (1/sig2u^2)*(Xstar %*% Hinv) %*% t(Xstar)
  }
  # MATRIX C IN (12). PAPER SAP
  if (env$type == "sar") {
    der2_reml_rho <- - t(y) %*%
      (t(Wsp) %*% (P %*% (Wsp  %*% y))) -
      sum(diag(solve(A1, Wsp)^2))
    der2_reml_rho <- as.numeric(der2_reml_rho )
    } else { der2_reml_rho <- NULL } 
  return(der2_reml_rho)
}
  
  


