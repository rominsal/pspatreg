fit_pspat <- function(env, con) {
  y <- env$y
  nfull <- length(y)
  if (!is.null(env$Xfull)) {
    X <- Matrix(env$Xfull)
  } else X <- matrix(0, nrow = length(y), ncol = 1)
  if (!is.null(env$Zfull)) {
    Z <- Matrix(env$Zfull)
  } else Z <- matrix(0, nrow = length(y), ncol = 1)
  Wsp <- env$Wsp
  sp1 <- env$sp1
  sp2 <- env$sp2
  time <- env$time
  nfull <- env$nfull
  nsp <- env$nsp
  nt <- env$nt
  Ifull <- Diagonal(nfull)
  Isp <- Diagonal(nsp)
  It <- Diagonal(nt)
  if (is.null(con$vary_init)) {
    la <-  var(as.numeric(y)) 
  } else  la <- con$vary_init 
  namesla <- c("sig2u")
  # Build vector and matrices for variance components in mixed model
  if (is.null(time)) 
    var_comp <- par_var_comp2d(la = la, env)
  else var_comp <- par_var_comp3d(la = la, env)
  # Number of parameters
  np <- var_comp$np
  np_eff <- var_comp$np_eff
  assign("np", np, envir = env)
  assign("np_eff", np_eff, envir = env)
  # Vector of parameters
  la <- var_comp$la
  if (env$nvarspt > 0) {
    if (is.null(time)) { 
      if (!env$psanova) {
        namestauspt <- paste("tausp", 1:2, sep = "")
      } else {
        namestauspt <- c("tauf1_main", "tauf2_main", 
                        "tauf12.1",  "tauf12.2")
      } 
    } else { # nt > 1
      if (!env$psanova) {
        namestauspt <- c(paste("tausp", 1:2, sep = ""),
                         "tautime")
      } else {
        namestauspt <- c("tauf1_main", 
                         "tauf2_main", 
                         "tauft_main",
                         "tauf12.1", "tauf12.2",
                         "tauf1t.1", "tauf1t.2",
                         "tauf2t.1", "tauf2t.2",
                         "tauf12t.1", "tauf12t.2", 
                         "tauf12t.3")  
      }
    }
    namesla <- c(namesla, namestauspt)
  }
  if (env$nvarnopar > 0) {
    namestaunopar <- paste("taunopar", 1:env$nvarnopar, 
                         sep = "")
    namesla <- c(namesla, namestaunopar)
  }
  names(la) <- namesla
  # Do not remove var_comp, it is used in the next loop...
  # 0 0
  # 0 I
  D <- Matrix(diag(c(rep(0, np_eff[1]), 
                             rep(1, sum(np_eff[-1])))))
  assign("D", D, envir = env)
  env$D <- D
  ## initialize rho,  delta and phi
  if (env$type %in% c("sar", "sdm", "sarar")) 
    rho <- con$rho_init  
  else  rho <- NULL 
  if (env$type %in% c("sem", "sdem", "sarar")) 
    delta <- con$delta_init  
  else  delta <- NULL 
  if (env$cor == "ar1" ) 
    phi <- con$phi_init  
  else  phi <- NULL 
  # Initialise the parameters
  if (is.null(con$bold)) bold = rep(0, sum(np_eff))
  if (length(np_eff) > 1) {
    eta <- X %*% bold[1:np_eff[1]] + 
      Z %*% bold[-(1:np_eff[1])] #+ offset
  } else eta <- X %*% bold[1:np_eff[1]] # Only fixed effects
  start <- proc.time()[3]
  for (iq in 1:con$maxit) {
    # Nested loops for spatial parameters and SAP
    if (!is.null(rho)) {
      A1 <- Isp - rho*Wsp } else if (nsp > 1) { 
        A1 <- Isp } else A1 <- Ifull 
    if (!is.null(delta)) {
      A2 <- Isp - delta*Wsp } else if (nsp > 1) { 
        A2 <- Isp } else A2 <- Ifull
    if (!is.null(time) || nt > 1) {
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
      # ystar_check <- matrix(kronecker(A2 %*% A1, 
      #                                 LcholOmegainv) %*% y, 
      #                       ncol = 1)
      #range(ystar - ystar_check)
      ## COMPROBAR SI Xstar Y Zstar PUEDEN CALCULARSE ASÍ...
      # Xstar <- apply(aperm(RH(Isp, RH(LcholOmegainv,
      #                               array(X, 
      #                    dim=c(ncol(LcholOmegainv), 
      #                          ncol(Isp), ncol(X))))),
      #                    perm=c(2, 1, 3)), 2, rbind)
      # Zstar <- apply(aperm(RH(Isp, RH(LcholOmegainv,
      #                              array(Z, 
      #                    dim=c(ncol(LcholOmegainv), 
      #                              ncol(Isp), ncol(Z))))),
      #                    perm=c(2, 1, 3)), 2, rbind)
      Xstar <- kronecker(A2, LcholOmegainv) %*% X
      Zstar <- kronecker(A2, LcholOmegainv) %*% Z
    } else { 
      Omega <- Omegainv <- LcholOmegainv <- It
      ystar <- A2 %*% (A1 %*% y)
      Xstar <- A2 %*% X
      Zstar <- A2 %*% Z
    }
    edftot_old <- 0
    for (it in 1:con$maxit) {
      if (length(np_eff) > 1) {
          # Model includes random effects
          # Builds covariance G matrix for random effects
          if (is.null(time)) {
            lG <- build_G2d(la = la, lg = var_comp, env)
          } else {
            lG <- build_G3d(la = la, lg = var_comp, env)
          }
          G <- lG$G
          Ginv <- lG$Ginv
          G_eff <- lG$G_eff
          Ginv_eff <- lG$Ginv_eff
          assign("G", G, envir = env)
          assign("Ginv", Ginv, envir = env)
          assign("G_eff", G_eff, envir = env)
          assign("Ginv_eff", Ginv_eff, envir = env)
          rm(lG)
          sig2u <- la["sig2u"]
          assign("sig2u", sig2u, envir = env)
          #if (is.null(weights)) {
          #    w <- as.vector(matrix(1,nrow=length(y))) } else { w <- weights }
          mat <- construct_matrices(Xstar, Zstar, ystar)
          # MATRIX C IN (12). PAPER SAP
          C <- Matrix( 
            construct_block(mat$XtX, 
                            t(mat$ZtX*G_eff), 
                            mat$ZtX, 
                            t(mat$ZtZ*G_eff)))
          H <- (1/sig2u)*C + D
          Hinv <- try(solve(H))
          if (inherits(Hinv, "try-error"))
            Hinv <- ginv(as.matrix(H))
          b <- as.vector((1/sig2u)*Hinv %*% mat$u)
          bfixed <- b[1:np_eff[1]]
          names(bfixed) <- gsub("X_", "", colnames(X))
          names(bfixed) <- paste("fixed_", 
                                 names(bfixed), sep = "")
          assign("bfixed", bfixed, envir = env)
          brandom <- G_eff*b[-(1:np_eff[1])]
          names(brandom) <- gsub("Z_", "", colnames(Z))
          names(brandom) <- paste("random_", 
                                  names(brandom), sep = "")
          assign("brandom", brandom, envir = env)
          # Compute effective dimensions and variances
          # Only the diagonal of ZtPZ
          dZtPZ <- 1/la[1] * apply((t(Hinv[-(1:np_eff[1]), ]) * 
                     mat$ZtXtZ), 2, sum)
          ## Check (8). Paper SAP
           # n1 <- ncol(env$Zfull)
           # n2 <- ncol(env$Xfull)
           # pr <- cbind(
           #      matrix(0, nrow = n1,
           #               ncol = n2),
           #      Diagonal(n1)
           #      )
           # dZtPZ2 <- (1/la[1]*pr %*%
           #   Hinv %*% mat$ZtXtZ)
           # diag_dZtPZ2 <- diag(as.matrix(dZtPZ2))
           # range(dZtPZ - diag_dZtPZ2)
          index.zeros.G <- G == 0
          brandom_wide <- dZtPZ_wide <- rep(0, length(G))
          brandom_wide[!index.zeros.G] <- brandom
          dZtPZ_wide[!index.zeros.G] <- dZtPZ
          if (is.null(time)) {
            ltau_edf <- update_tau2d(la = la, 
                                     lg = var_comp, 
                                     G = G,
                                     dZtPZ_wide = dZtPZ_wide, 
                                     brandom_wide = brandom_wide,
                                     env = env)
          } else {
            ltau_edf <- update_tau3d(la = la, 
                                     lg = var_comp, 
                                     G = G,
                                     dZtPZ_wide = dZtPZ_wide, 
                                     brandom_wide = brandom_wide,
                                     env = env)
          }
          if (env$nvarspt > 0) { 
            if (is.null(time)) {
              tau1 <- ltau_edf$tau1
              tau2 <- ltau_edf$tau2
              ed1 <- ltau_edf$ed1
              ed2 <- ltau_edf$ed2
              tauspt <- c(tau1, tau2)
              edfspt <- c(ed1, ed2)
              if (env$psanova) {
                tau3 <- ltau_edf$tau4
                tau4 <- ltau_edf$tau4
                ed3 <- ltau_edf$ed3
                ed4 <- ltau_edf$ed4
                tauspt <- c(tauspt, tau3, tau4)
                edfspt <- c(edfspt, ed3, ed4)
                }
              } else { 
                tau1 <- ltau_edf$tau1
                tau2 <- ltau_edf$tau2
                tau3 <- ltau_edf$tau3
                ed1 <- ltau_edf$ed1
                ed2 <- ltau_edf$ed2
                ed3 <- ltau_edf$ed3
                tauspt <- c(tau1, tau2, tau3)
                edfspt <- c(ed1, ed2, ed3)
                if (env$psanova) {
                  tau4 <- ltau_edf$tau4 
                  tau5 <- ltau_edf$tau5
                  tau6 <- ltau_edf$tau6
                  tau7 <- ltau_edf$tau7 
                  tau8 <- ltau_edf$tau8 
                  tau9 <- ltau_edf$tau9
                  tau10 <- ltau_edf$tau10
                  tau11 <- ltau_edf$tau11 
                  tau12 <- ltau_edf$tau12
                  ed4 <- ltau_edf$ed4
                  ed5 <- ltau_edf$ed5
                  ed6 <- ltau_edf$ed6
                  ed7 <- ltau_edf$ed7
                  ed8 <- ltau_edf$ed8
                  ed9 <- ltau_edf$ed9
                  ed10 <- ltau_edf$ed10
                  ed11 <- ltau_edf$ed11
                  ed12 <- ltau_edf$ed12
                  tauspt <- c(tauspt, tau4, tau5, tau6, tau7, 
                              tau8, tau9, tau10, tau11, tau12)
                  edfspt <- c(edfspt, ed4, ed5, ed6, ed7, ed8, 
                              ed9, ed10, ed11, ed12)
                }
              }
            names(tauspt) <- namestauspt
            names(edfspt) <- gsub("tau", "", namestauspt)
          } else {
            tauspt <- NULL
            edfspt <- NULL
          }
          if (env$nvarnopar > 0) {
            taunopar <- ltau_edf$taunopar
            # Add 1 for fixed effects of each nonparametric variable
            # edfnopar <- ltau_edf$edfnopar + 1
            edfnopar <- ltau_edf$edfnopar
            names(edfnopar) <- paste("edfnopar", 1:env$nvarnopar, 
                                     sep = "")
          } else {
            taunopar <- NULL
            edfnopar <- NULL
          }
          # Regression (Fahrmeir et al.) pp. 180
          ssr <- as.numeric(mat$yty - t(c(bfixed, brandom)) %*% 
                              (2*mat$u - C %*% b))
          ## Checking ssr
          # resids <- y - Xfull %*% bfixed - Zfull %*% brandom
          # ssr2 <- sum(resids^2)
          # New variance
          lanew <- c(sig2u)
          namesX <- colnames(X)
          namesX.nopar <- namesX[grepl("pspl", namesX)]
          # Fixed effects of X.nopar has been added to edfnopar.
          edftot <- length(namesX) - length(namesX.nopar)
          if (env$nvarspt > 0) {
            lanew <- c(lanew, tauspt)
            edftot <- edftot + sum(edfspt) 
          }
          if (env$nvarnopar > 0) {
            lanew <- c(lanew, taunopar)
            edftot <- edftot + sum(edfnopar)
          }
        } else { # Model without random effects
          tauspt <- taunopar <- NULL
          edfspt <- edfnopar <- NULL
          sig2u <- la["sig2u"]
          assign("sig2u", sig2u, envir = env)
          mat <- construct_matrices(Xstar, Zstar, ystar)
          # MATRIX C IN (12). PAPER SAP
          C <- Matrix(mat$XtX)
          H <- (1/sig2u)*C 
          Hinv <- try(solve(H))
          if (inherits(Hinv, "try-error"))
            Hinv <- ginv(as.matrix(H))
          b <- as.vector((1/sig2u)*Hinv %*% mat$u[1:np_eff[1]])
          bfixed <- b[1:np_eff[1]]
          names(bfixed) <- gsub("X_", "", colnames(X))
          names(bfixed) <- paste("fixed_", 
                                 names(bfixed), sep = "")
          assign("bfixed", bfixed, envir = env)
          brandom <- 0
          assign("brandom", brandom, envir = env)
          ssr <- as.numeric(mat$yty - t(c(bfixed)) %*% 
                              (2*mat$u[1:np_eff[1]] - C %*% b))
          lanew <- c(sig2u)
          edftot <- ncol(X)
        } # end  Model without random effects
	      if (!is.null(rho)) edftot <- edftot + 1 
	      if (!is.null(delta)) edftot <- edftot + 1 
	      sig2u <- as.numeric((ssr/(length(y) - edftot)))
	      # Update first component of la with new sig2u
	      lanew["sig2u"] <- sig2u
	      assign("sig2u", sig2u, envir = env)
	      # Update rho and delta
	      dla <- mean(abs(la - lanew))
	      la <- lanew
	      dedftot <- edftot - edftot_old
	      edftot_old <- edftot
	      if (con$trace) {
	        message(paste('\n Iteration SAP: ', it))
	        message(paste('\n mean change in penalization param.: ', dla))
	        message(paste('\n change in edftot: ', dedftot))
	        message(paste('\n sig2u ', la[1]))
	        if (env$nvarspt > 0)
	          message(paste('\n edfspt:', round(edfspt, 2)))
	        if (env$nvarnopar > 0) 
	          message(paste('\n edfnopar: ', 
	              round(edfnopar, 2))) }
	      #  convergence check
	      if (nfull > 500) {
	        if (dedftot < con$tol2) break
	      } else {
	        if (dla < con$tol1) break
	      } 
	  } # end for (it in 1:maxit)
    if (!(env$type %in% c("sim", "slx")) || env$cor == "ar1") { 
      # Optimize using spatial and/or correlation parameters
      param <- namesparam <- NULL
      lower_par <- upper_par <- NULL
      if (!is.null(rho)) {
        param <- c(param, rho)
        namesparam <- c(namesparam, "rho") 
        lower_par <- c(lower_par, env$interval[1])
        upper_par <- c(upper_par, env$interval[2])
      } 
      if (!is.null(delta)) {
        param <- c(param, delta)
        namesparam <- c(namesparam, "delta")
        lower_par <- c(lower_par, env$interval[1])
        upper_par <- c(upper_par, env$interval[2])
      } 
      if (!is.null(phi)) {
        param <- c(param, phi)
        namesparam <- c(namesparam, "phi")
        lower_par <- c(lower_par, -0.99)
        upper_par <- c(upper_par, 0.99)
      } 
      names(param) <- namesparam
      if (con$optim == "llik_reml") {
        par_optim <- bobyqa(par = param, 
                            fn = llikc_reml,
                            lower = lower_par,
                            upper = upper_par,
                            control = list(rhobeg = 0.5,
                                           iprint = 0),
                            env = env)
        param_new <- par_optim$par
      }
      if (con$optim == "llik") {
        par_optim <- bobyqa(par = param, 
                            fn = llikc,
                            lower = lower_par,
                            upper = upper_par,
                            control = list(rhobeg = 0.5,
                                           iprint = 0),
                            env = env)
        param_new <- par_optim$par 
      }
      names(param_new) <- namesparam
      dparam <- mean(abs(param - param_new))
      param <- param_new
      if (!is.null(rho)) rho <- param["rho"]
      if (!is.null(delta)) delta <- param["delta"]
      if (!is.null(phi)) phi <- param["phi"]
    } else {  
      rho <- delta <- phi <- NULL
      param <- c(rho, delta, phi)
      param_optim <- param
      dparam <- 0
    }
    if (con$trace) {
      message(paste("\n Iteration Spatials and/or 
                    Correlation Parameters: ", iq))
      if (!is.null(rho)) 
        message(paste("\n  rho: ", rho))
      if (!is.null(delta)) 
        message(paste("\n  delta: ", delta))
      if (!is.null(phi)) 
        message(paste("\n  phi: ", phi))
    }
    # Check Convergence 
	  if (dparam < con$tol3) break
	} # End loop 
  end <- proc.time()[3]
	message(paste("\n Time to fit the model: ", round(end - start, 2), 
	    "seconds \n"))
	param_optim <- param
	if (length(np_eff) > 1) {
	  eta <- X %*% bfixed + Z %*% brandom #+ offset
	} else eta <- X %*% bfixed # Only fixed effects
#  FINAL ESTIMATES OF PARAMETERS
  #sig2u <- la["sig2u"]
  #assign("sig2u", sig2u, envir = env)
  assign("tauspt", tauspt, envir = env)
  assign("taunopar", taunopar, envir = env)
  assign("edfspt", edfspt, envir = env)
  assign("edfnopar", edfnopar, envir = env)
  assign("edftot", edftot, envir = env)
  # Valor de log.lik y log.lik.reml en el óptimo
  llikc_reml_optim <- -llikc_reml(param_optim, env)
  llikc_optim <- -llikc(param_optim, env)
  if (!(env$type %in% c("sim", "slx")) || env$cor == "ar1") {
    hessian_optim <- hessian(llikc_reml, 
                             param_optim, env = env)
    var_num <- solve(hessian_optim) 
    se_num <- sqrt(diag(var_num))
    names(var_num) <- names(se_num) <- namesparam
  } else var_num <- se_num <- NULL
  # se_an <- NULL
  # CHANGE WHEN IT IS READY anhess_llikc_reml FUNCTION...
  # if (env$type == "sar" && !(con$fdHess)) { 
  #   se_an <- sqrt(-(1/anhess_llikc_reml_2d(param_optim, 
  #                                              env)))
  #   names(se_an) <- namesparam
  #   } else  se_an <- NULL
  se_rho <- se_num["rho"]
  se_delta <- se_num["delta"]
  se_phi <- se_num["phi"]
  
  ########## COVARIANCE MATRICES FIXED AND RANDOM EFFECTS
  ## pp.375 Fahrmeir et al.
  ## Bayesian Covariance Matrix
  if (con$trace) start <- proc.time()[3]
  if (!is.null(rho)) {
    A1 <- Isp - rho*Wsp } else { A1 <- Isp } 
  if (!is.null(delta)) {
    A2 <- Isp - delta*Wsp } else { A2 <- Isp }
  if (!is.null(time) || nt > 1) {
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
  } else { 
    Omega <- Omegainv <- LcholOmegainv <- It
    ystar <- Matrix(A2 %*% (A1 %*% y))
    Xstar <- Matrix(A2 %*% X)
    Zstar <- Matrix(A2 %*% Z)
  }
  XZstar <- cbind(Xstar, Zstar)
  if (length(np_eff) > 1) {
    DG <- Matrix(diag(c(rep(0, np_eff[1]), Ginv_eff)))
  } else { # Only random effects
    XZstar <- XZstar[, 1:np_eff]
    DG <- 0
  }
  XZt_Rinv_XZ <- (1/sig2u)*crossprod(XZstar)
  Var_By <- as(solve(XZt_Rinv_XZ + DG, tol = 1e-30),
               "dpoMatrix")
  rownames(Var_By) <- colnames(Var_By) <- c(names(bfixed), 
                                            names(brandom))
  seby_bfixed <- sqrt(diag(as.matrix(Var_By[names(bfixed), 
                                            names(bfixed)])))
  names(seby_bfixed) <- names(bfixed)
  if (length(np_eff) > 1) {
    seby_brandom <- sqrt(diag(as.matrix(Var_By[names(brandom),
                                               names(brandom)])))
    names(seby_brandom) <- names(brandom)
  } else seby_brandom <- NULL
  ## Frequentist Covariance Matrix
  Var_Fr <- Var_By %*% XZt_Rinv_XZ %*% Var_By
  rownames(Var_Fr) <- colnames(Var_Fr) <- c(names(bfixed), 
                                            names(brandom))
  sefr_bfixed <- sqrt(diag(as.matrix(Var_Fr[names(bfixed), 
                                            names(bfixed)])))
  names(sefr_bfixed) <- names(bfixed)
  if (length(np_eff) > 1) {
    sefr_brandom <- sqrt(diag(as.matrix(Var_Fr[names(brandom),
                                               names(brandom)])))
    names(sefr_brandom) <- names(brandom)
  } else sefr_brandom <- NULL
  
  if (con$trace) {
    end <- proc.time()[3]
    message(paste("\n Time to compute covariances: ", 
        (end-start), "seconds")) }
  # Fits and Resids
  if (length(np_eff) > 1) {
    fit_A1y <- X %*% bfixed + Z %*% brandom # + offset
    XZ <- cbind(X, Z)
  } else {
    fit_A1y <- X %*% bfixed
    XZ <- X # Only fixed effects
  }
  if (is.null(time) && nt == 1) {
    fit <- solve(A1, fit_A1y)
    seby_fit_A1y <- rowSums((XZ %*% Var_By) * XZ)^0.5
    seby_fit <- rowSums((solve(A1, XZ) %*% Var_By) *
                                  solve(A1, XZ))^0.5
    sefr_fit_A1y <- rowSums((XZ %*% Var_Fr) * XZ)^0.5
    sefr_fit <- rowSums((solve(A1, XZ) %*% Var_Fr) *
                                  solve(A1, XZ))^0.5
    #residuals <- as.vector((A1 %*% y) - fit_A1y)
  } else {  
    fit <- solve( kronecker(A1, It), fit_A1y)
    seby_fit_A1y <- rowSums((XZ %*% Var_By) * XZ)^0.5
    seby_fit <- rowSums((
      solve( kronecker(A1, It), XZ) %*% Var_By) *
      solve( kronecker(A1, It), XZ))^0.5
    sefr_fit_A1y <- rowSums((XZ %*% Var_Fr) * XZ)^0.5
    sefr_fit <- rowSums((
      solve( kronecker(A1, It), XZ) %*% Var_Fr) *
      solve( kronecker(A1, It), XZ))^0.5
    # residuals <- as.vector((
    #   kronecker(A1, It) %*% y) - fit_A1y)
  }
  residuals <- as.vector(y - fit)
  # Compute AIC y BIC based on loglik functions 
  # (Fahrmeir, pp. 664 and 677)
  aic <- -2*llikc_optim + 2*edftot
  bic <- -2*llikc_optim + log(length(y))*edftot
  res <- list(edfspt = edfspt, 
              edfnopar = edfnopar, 
              edftot = edftot,
              tauspt = tauspt, 
              taunopar = taunopar,
              fitted.values = as.vector(fit),
              fit_A1y = as.vector(fit_A1y),
              seby_fitted.values = as.vector(seby_fit),
              seby_fit_A1y = as.vector(seby_fit_A1y),
              sefr_fitted.values = as.vector(sefr_fit),
              sefr_fit_A1y = as.vector(sefr_fit_A1y),
              residuals = as.vector(residuals),
              sig2 = sig2u,
              rho = rho,
              se_rho = se_rho,
              delta = delta,
              se_delta = se_delta,
              phi = phi,
              se_phi = se_phi,
              bfixed = bfixed, 
              brandom = brandom,
              seby_bfixed = seby_bfixed, 
              seby_brandom = seby_brandom,
              sefr_bfixed = sefr_bfixed, 
              sefr_brandom = sefr_brandom,
              llik = llikc_optim, 
              llik_reml = llikc_reml_optim,
              aic = aic, 
              bic = bic,
              vcov_by = Var_By,
              vcov_fr = Var_Fr,
              sp1 = env$sp1, 
              sp2 = env$sp2, 
              time = env$time,
              psanova = env$psanova)

} # end of function

